#ifndef VIENNASHE_SHE_ASSEMBLE_SCATTERING_HPP
#define VIENNASHE_SHE_ASSEMBLE_SCATTERING_HPP

/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include <cmath>

// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/assemble_common.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/util/filter.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"

#include "scattering/all.hpp"

/** @file viennashe/she/assemble_scattering.hpp
    @brief Generic assembly of the scattering operator(s) is implemented here
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Finds the closest H-index for the provided kinetic energy
     *
     * @param quan             The SHE quantities
     * @param elem                   The vertex/edge on which to search
     * @param index_H_hint           Index of the energy from which to scatter
     * @param kinetic_energy         Kinetic energy for which index_H should be found
     *
     * @return Indices in vector 'energies' that are closest to the lower and upper energy after scattering.
     */
    template <typename SHEQuantity, typename ElementType>
    long energy_index_H(SHEQuantity const & quan,
                        ElementType const & elem,
                        long index_H_hint,
                        double kinetic_energy)
    {
      long result_index = index_H_hint;

      if (   kinetic_energy < quan.get_kinetic_energy(elem, std::size_t(result_index))
          || kinetic_energy > quan.get_kinetic_energy(elem, std::size_t(result_index)) ) // this is an inelastic process, so we have to search
      {
        long increment = 1;

        double diff = kinetic_energy - quan.get_kinetic_energy(elem, std::size_t(result_index));
        double diff_init = diff;

        if (diff < 0) //we need to search for a lower index
          increment = -1;

        while ( (result_index < static_cast<long>(quan.get_value_H_size() - 1)) && (result_index > 0)
               && (diff * diff_init > 0) )  //sign hasn't changed yet
        {
          result_index += increment;
          diff = kinetic_energy - quan.get_kinetic_energy(elem, std::size_t(result_index));
        }

        if (diff * diff_init > 0) //still no sign change, hence final scattering state is outside the considered range
          return -1;

        // Pick best approximation to final energy:
        double error_1 = std::abs(diff);
        double error_2 = error_1;
        if ( (result_index - increment >= 0)
            && (result_index - increment < static_cast<long>(quan.get_value_H_size())) )
          error_2 = std::abs(kinetic_energy - quan.get_kinetic_energy(elem, std::size_t(result_index - increment)));

        if (error_2 < error_1)
          return result_index - increment;  //going one step/increment back gives the best approximation
      }

      return result_index;
    }


    template <typename ScatterProcessesT,
              typename DeviceType,
              //typename TimeStepQuantitiesT,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType,
              typename ElementType,
              typename CouplingMatrix>
    void assemble_scattering_operator_on_box(ScatterProcessesT const & scatter_processes,
                                             DeviceType const & device,
                                             viennashe::config const & conf,
                                             SHEQuantity const & quan,
                                             MatrixType & A, VectorType & b,
                                             ElementType const & elem, std::size_t index_H,
                                             CouplingMatrix const & coupling_in_scatter,
                                             CouplingMatrix const & coupling_out_scatter)
    {
      typedef typename viennagrid::result_of::point<typename DeviceType::mesh_type>::type   PointType;

      typedef scattering_base<DeviceType>    ScatterProcessType;

      //
      // Step 1: Get standard quantities
      //
      const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

      const double pi = viennashe::math::constants::pi;

      const long row_index = quan.get_unknown_index(elem, index_H);

      if (row_index < 0)
        return;

      double kinetic_energy = quan.get_kinetic_energy(elem, index_H);

      typedef typename ScatterProcessType::value_type    scatter_descriptors_type;

      const long expansion_order_mid = static_cast<long>(quan.get_expansion_order(elem, index_H));

      bool odd_assembly = detail::is_odd_assembly(elem, viennagrid::cells(device.mesh())[0]);

      double box_volume = viennagrid::volume(elem);
      if (odd_assembly)
        box_volume *= detail::cell_connection_length(device.mesh(), elem, viennagrid::cells(device.mesh())[0]) / static_cast<double>(PointType::dim);

      //
      // Now iterate over all scattering operators
      //
      for (typename ScatterProcessesT::const_iterator scatter_process_it  = scatter_processes.begin();
                                                      scatter_process_it != scatter_processes.end();
                                                    ++scatter_process_it)
      {
        ScatterProcessType const & scatter_process = *(*scatter_process_it);

        if (log_assemble_scattering_operator::enabled)
          log::debug<log_assemble_scattering_operator>() << "* assemble_scattering_operator_on_box(): Assembling scattering process with ID " << scatter_process.id() << std::endl;

        scatter_descriptors_type scatter_events = scatter_process(elem, kinetic_energy, quan.get_carrier_type_id());

        //
        // Part 2: Loop over all possible scattering events for this energy (in- and out-scattering)
        //
        for (std::size_t scatter_index = 0; scatter_index < scatter_events.size(); ++scatter_index)
        {
          const double initial_kinetic_energy = scatter_events[scatter_index].initial_energy();
          const double   final_kinetic_energy = scatter_events[scatter_index].final_energy();
          const double rate = scatter_events[scatter_index].rate();
          const double generation_rate = scatter_events[scatter_index].generation_rate();

          // skip if no contribution
          if (rate <= 0 && generation_rate <= 0)  // inequality in order to avoid strict comparison with zero (triggers warnings)
            continue;

          long initial_state_index_H = energy_index_H(quan, elem, static_cast<long>(index_H), initial_kinetic_energy);
          long   final_state_index_H = energy_index_H(quan, elem, static_cast<long>(index_H),   final_kinetic_energy);


          // skip this process if one of the states is outside the simulation region
          if ( (initial_state_index_H == -1) || (final_state_index_H == -1))
            continue;

          //
          // Compute density of states at initial and final states
          //
          const double Z_initial = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), elem, static_cast<std::size_t>(initial_state_index_H));
          const double Z_final   = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), elem, static_cast<std::size_t>(final_state_index_H));

          const double energy_height = std::min(box_height(quan, elem, static_cast<std::size_t>(initial_state_index_H)),  // Note: box height near band edge might be smaller, reducing the active area.
                                                box_height(quan, elem, static_cast<std::size_t>(final_state_index_H)));

          const bool elastic_scattering_roundoff_error_prevention = (initial_state_index_H == final_state_index_H) && !odd_assembly;

          //
          // Out-Scattering
          //
          if (initial_kinetic_energy <= kinetic_energy && initial_kinetic_energy >= kinetic_energy)
          {
            const long col_index = quan.get_unknown_index(elem, static_cast<std::size_t>(final_state_index_H));

            if (col_index >= 0) //only write contribution if final state is within simulation domain
            {
              viennashe::math::harmonics_iteration_type harmonics_it_id = (odd_assembly ? viennashe::math::ODD_HARMONICS_ITERATION_ID : viennashe::math::EVEN_HARMONICS_ITERATION_ID);
              double coefficient = 0;
              switch (conf.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                coefficient = 2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height;
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                coefficient = 2.0 * pi * rate * Z_final * box_volume * energy_height;
                break;
              default:
                throw std::runtime_error("assemble_scattering_operator_on_box(): Unknown SHE discretization type!");
              }
              viennashe::util::add_block_matrix(A,
                                                static_cast<std::size_t>(row_index), static_cast<std::size_t>(row_index),
                                                coefficient,
                                                coupling_out_scatter,
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                elastic_scattering_roundoff_error_prevention
                                              );

              if (log_assemble_scattering_operator::enabled)
              {
                log::debug<log_assemble_scattering_operator>() << "* assemble_scattering_operator_on_box(): out-scattering contribution: " << std::endl;
                log::debug<log_assemble_scattering_operator>() << "      2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height " << std::endl;
                log::debug<log_assemble_scattering_operator>() << "   at (" << row_index << ", " << row_index << ")" << std::endl;
                log::debug<log_assemble_scattering_operator>() << "            rate = " << rate << std::endl;
                log::debug<log_assemble_scattering_operator>() << "       Z_initial = " << Z_initial << std::endl;
                log::debug<log_assemble_scattering_operator>() << "         Z_final = " << Z_final << std::endl;
                log::debug<log_assemble_scattering_operator>() << "      box_volume = " << box_volume << std::endl;
                log::debug<log_assemble_scattering_operator>() << "   energy_height = " << energy_height << std::endl;
              }

              if (with_full_newton)
              {
                // We assume 'rate' to be independent of the potential, which is certainly true for phonon scattering, but also for most other processes.
                // Also note that dQ/dphi is zero, hence no additional entries in the Jacobian matrix.

                // residual contribution:
                throw std::runtime_error("assemble_scattering_operator_on_box(): Newton not available!");
                viennashe::util::subtract_folded_block_vector(b,
                                                              static_cast<std::size_t>(row_index),
                                                              2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height,
                                                              coupling_out_scatter,
                                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                              quan.get_values(elem, index_H));
              }

            }
          }


          // Note: A scattering process with initial_kinetic_energy == final_kinetic_energy leads to both in- and out-scattering. Thus, do not try to merge the two if-statements above and below!

          //
          // In-Scattering
          //
          if (final_kinetic_energy <= kinetic_energy && final_kinetic_energy >= kinetic_energy) //in-scattering
          {
            const long col_index = quan.get_unknown_index(elem, static_cast<std::size_t>(initial_state_index_H));

            if (col_index >= 0)  //If final scatter state is invalid, no scattering considered!
            {
              viennashe::math::harmonics_iteration_type harmonics_it_id = (odd_assembly ? viennashe::math::ODD_HARMONICS_ITERATION_ID : viennashe::math::EVEN_HARMONICS_ITERATION_ID);
              double coefficient = 0;
              switch (conf.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                coefficient = -2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height
                              - generation_rate * Z_initial * box_volume * energy_height;  // additional generation due to in-scattering
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                coefficient = -2.0 * pi * rate * Z_final * box_volume * energy_height
                              - generation_rate * box_volume * energy_height;  // additional generation due to in-scattering
                break;
              default:
                throw std::runtime_error("assemble_scattering_operator_on_box(): Unknown SHE discretization type!");
              }
              viennashe::util::add_block_matrix(A,
                                                static_cast<std::size_t>(row_index), static_cast<std::size_t>(col_index),
                                                coefficient,                           // in-scattering
                                                coupling_in_scatter,
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                elastic_scattering_roundoff_error_prevention
                                               );

              if (log_assemble_scattering_operator::enabled)
              {
                log::debug<log_assemble_scattering_operator>() << "* assemble_scattering_operator_on_box(): in-scattering contribution:" << std::endl;
                log::debug<log_assemble_scattering_operator>() << "      - 2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height" << std::endl;
                log::debug<log_assemble_scattering_operator>() << "      - generation_rate * Z_initial * box_volume * energy_height" << std::endl;
                log::debug<log_assemble_scattering_operator>() << "  at (" << row_index << ", " << col_index << ")" << std::endl;
                log::debug<log_assemble_scattering_operator>() << "              rate = " << rate << std::endl;
                log::debug<log_assemble_scattering_operator>() << "         Z_initial = " << Z_initial << std::endl;
                log::debug<log_assemble_scattering_operator>() << "           Z_final = " << Z_final << std::endl;
                log::debug<log_assemble_scattering_operator>() << "        box_volume = " << box_volume << std::endl;
                log::debug<log_assemble_scattering_operator>() << "     energy_height = " << energy_height << std::endl;
                log::debug<log_assemble_scattering_operator>() << "   generation_rate = " << generation_rate << std::endl;
              }

              if (with_full_newton)
              {
                // We assume 'rate' to be independent of the potential, which is certainly true for phonon scattering, but also for most other processes.
                // Also note that dQ/dphi is zero, hence no additional entries in the Jacobian matrix.

                // residual contribution:
                throw std::runtime_error("assemble_scattering_operator_on_box(): Newton not available!");
                viennashe::util::subtract_folded_block_vector(b,
                                                              static_cast<std::size_t>(row_index),
                                                              - 2.0 * pi * rate * Z_initial * Z_final * box_volume * energy_height
                                                                - generation_rate * Z_initial * box_volume * energy_height,  // additional generation due to in-scattering
                                                              coupling_in_scatter,
                                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, harmonics_it_id),
                                                              quan.get_values(elem, index_H));
              } // with_full_newton

            } // if col_index >= 0
          } // in-scattering

        } // for all scattering events (i.e. events within a process)

      } //for all scattering processes
    }


  } //namespace she
} //namespace viennashe

#endif
