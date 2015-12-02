#ifndef VIENNASHE_SHE_ASSEMBLE_TIMEDERIVATIVE_HPP
#define VIENNASHE_SHE_ASSEMBLE_TIMEDERIVATIVE_HPP


/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include "viennashe/forwards.h"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/she/assemble_common.hpp"
#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/filter.hpp"

namespace viennashe
{
  namespace she
  {
    namespace detail
    {
      /** @brief Obtain the potential difference from an edge */
      template <typename ElementType, typename QuantityPotential, typename QuantityPotentialOld>
      double potential_difference_to_old_timestep(ElementType          const & elem,
                                                  QuantityPotential    const & new_potential,
                                                  QuantityPotentialOld const & old_potential)
      {
        return 0.0;//   0.5*(new_potential.get_value(viennagrid::vertices(elem)[0]) + new_potential.get_value(viennagrid::vertices(elem)[1]))
               //- 0.5*(old_potential.get_value(viennagrid::vertices(elem)[0]) + old_potential.get_value(viennagrid::vertices(elem)[1]));
      }

      /** @brief Obtain the potential difference from a vertex */
      template <typename ConfigType, typename QuantityPotential, typename QuantityPotentialOld>
      double potential_difference_to_old_timestep(viennagrid_element_id        elem,
                                                  QuantityPotential    const & new_potential,
                                                  QuantityPotentialOld const & old_potential)
      {
        return new_potential.get_value(elem) - old_potential.get_value(elem);
      }
    }


    template <typename DeviceType,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType,
              typename ElementType,
              typename CouplingMatrixType,
              typename QuantityPotential,
              typename QuantityPotentialOld
              >
    void assemble_timederivative( DeviceType const & device,
                                  viennashe::config const & conf,
                                  SHEQuantity const & quan,
                                  SHEQuantity const & quan_old,
                                  MatrixType & A,
                                  VectorType & b,
                                  ElementType const & elem, std::size_t index_H,
                                  CouplingMatrixType const & identity_matrix,
                                  QuantityPotential const & quan_pot,
                                  QuantityPotentialOld const & quan_pot_old,
                                  bool odd_assembly)
    {
      (void)quan_pot; (void)quan_pot_old; //avoid unused parameter warnings

      const double rtime = (conf.time_step_size() > 0.0) ? 1.0 / conf.time_step_size() : 0.0;
      const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);


      if (rtime == 0.0)
        return; //stationary

      const long row_index = quan.get_unknown_index(elem, index_H);
      if (row_index < 0)
        return; //no unknown here

      //if (odd_assembly)
      //  return;

      viennashe::math::harmonics_iteration_type harmonics_it_id = (odd_assembly ? viennashe::math::ODD_HARMONICS_ITERATION_ID : viennashe::math::EVEN_HARMONICS_ITERATION_ID);

      if (! with_full_newton )
      {
        const long expansion_order = static_cast<long>(quan.get_expansion_order(elem, index_H));
        double box_volume;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_volume(device.mesh(), elem, &box_volume));
        //if (odd_assembly)
        //  box_volume *= detail::cell_connection_length(device.mesh(), elem, viennagrid::cells(device.mesh())[0]) / PointType::dim;

        //const double Z               = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), elem, index_H);

        //
        // Assemble the df/dH term: \pm q * dphi/dt * f(H_{n \pm 1}, t_{k+1})
        //
        /*
        const double dphi = detail::potential_difference_to_old_timestep(elem, quan_pot, quan_pot_old); // new - old
        const double q    = viennashe::physics::constants::q;

        const viennashe::carrier_type_id carrier_id          = quan.get_carrier_type_id();
        const double                     polarity            = (carrier_id == ELECTRON_TYPE_ID) ? -1.0 : 1.0;

        const long col_index_minus = quan.get_unknown_index(elem, index_H - 1);
        if (col_index_minus >= 0 )
        {
          double Z_lower   = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), elem, index_H - 1);

          // from box integration
          const double coeff = polarity * Z_lower * box_volume;

          viennashe::util::add_block_matrix(A,
                                            row_index, col_index_minus,
                                            - 0.5 * q * coeff * dphi * rtime ,
                                            identity_matrix,
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                                          );
        }

        const long col_index_plus = quan.get_unknown_index(elem, index_H + 1);
        if (col_index_plus >= 0 )
        {
          double Z_upper   = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), elem, index_H + 1);

          const double coeff = polarity * Z_upper * box_volume;

          viennashe::util::add_block_matrix(A,
                                            row_index, col_index_plus,
                                            + 0.5 * q * coeff * dphi * rtime ,
                                            identity_matrix,
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                                          );
        }
        */

        //
        // df/dt term
        //
        const double height = box_height(quan, elem, index_H);
        const double coeff  = height * box_volume;

        viennashe::util::add_block_matrix(A,
                                          std::size_t(row_index), std::size_t(row_index),
                                          coeff * rtime,
                                          identity_matrix,
                                          viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                          viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                                        );

        viennashe::util::subtract_folded_block_vector(b,
                                                      std::size_t(row_index),
                                                      -coeff * rtime,
                                                      identity_matrix,
                                                      viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                                      viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                                      quan_old.get_values(elem,index_H));
      }
      else if (with_full_newton)
      {
        log::warning() << "* simulator(): !! WARNING !!  Time dependence is not available for Newton iterations! Skipping df/dt " << std::endl;
      }

    }



  } // she
} // viennashe

#endif /* ASSEMBLE_TIMEDERIVATIVE_HPP */

