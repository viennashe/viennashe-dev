#ifndef VIENNASHE_SHE_ASSEMBLE_TRAPS_HPP
#define VIENNASHE_SHE_ASSEMBLE_TRAPS_HPP

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


// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/setters.hpp"
#include "viennashe/accessors.hpp"
#include "viennashe/math/constants.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/assemble_common.hpp"
#include "viennashe/she/timestep_quantities.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/math/linalg_util.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"


#include "viennashe/models/srh_kinetics.hpp"


/** @file viennashe/she/assemble_traps.hpp
    @brief Generic assembly of traps is implemented here. At the moment only SRH type traps are supported
*/

namespace viennashe
{
  namespace she
  {

    namespace detail
    {

      /** @brief Assembles the coupling of the distribution function with the trap occupancy for a cell (even-order f_lm) */
      template <typename DeviceType,
                typename TimeStepQuantitiesT,
                typename SHEQuantity,
                typename MatrixType,
                typename VectorType,
                typename CouplingMatrixType>
      void assemble_traps_coupling_on_cell(DeviceType const & device,
                                              TimeStepQuantitiesT const & quantities,
                                              SHEQuantity const & quan,
                                              viennashe::config const & conf,
                                              typename DeviceType::cell_type const & el,
                                              std::size_t index_H,
                                              MatrixType & matrix, VectorType & rhs,
                                              CouplingMatrixType const & diagonal_coupling_matrix,
                                              CouplingMatrixType const & coupling_matrix_00
                                             )
      {
        typedef typename DeviceType::trap_level_container_type          TrapContainerType;
        typedef typename TrapContainerType::const_iterator              TrapIterator;

        const long row_index = quan.get_unknown_index(el, index_H);

        if (row_index < 0)
          return;

        TrapContainerType const & traps = device.get_trap_levels(el);

        const std::size_t num_trap_unknowns = quantities.num_trap_unknown_indices(el);

        if (num_trap_unknowns <= 0)
          return;

        if (num_trap_unknowns != traps.size())
          throw viennashe::invalid_value_exception("The number of traps configured in the device does not match the number of unknowns for traps!", static_cast<double>(num_trap_unknowns));

        const double Z                  = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), el, index_H);
        const long   expansion_order    = static_cast<long>(quan.get_expansion_order(el, index_H));
        const double energy_height      = box_height(quan, el, index_H);

        viennashe::math::harmonics_iteration_type harmonics_it_id = viennashe::math::EVEN_HARMONICS_ITERATION_ID;

        double volume_contribution = viennagrid::volume(el);

        std::size_t inner_index = 0;
        for ( TrapIterator trap_it = traps.begin();
                           trap_it != traps.end();
                         ++trap_it, ++inner_index)
        {
          log::debug<log_traps>() << el << " => " << *trap_it << std::endl;

          const double occupancy = quantities.trap_occupancy(el, inner_index);

          //
          // Term      + Gamma_rec N_T f                       (electrons)
          //             Gamma_rec N_T f                       (holes)
          //
          // Note: Gamma_rec should already include the trap occupancies !!!
          //

          const double trap_density = trap_it->density();
          // Note: Gamma_rec should already include the trap occupancies !!!
          const double gamma_recombination = viennashe::models::srh::gamma_recombination(*trap_it, device, el, conf,
                                                                                          quantities, quan.get_carrier_type_id(),
                                                                                          occupancy, index_H);

          viennashe::util::add_block_matrix(matrix,
                                            std::size_t(row_index), std::size_t(row_index),
                                            gamma_recombination * trap_density * Z * volume_contribution * energy_height,
                                            diagonal_coupling_matrix,
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                                          );

          //
          // Term      - (1 - f) Gamma_gen N_T       (electrons)
          //           - (1 - f) Gamma_gen N_T       (holes)
          //
          // Neglecting Pauli principle, thus only      - Gamma_gen N_T f_T       is assembled
          //
          // Note: Gamma_gen should already include the trap occupancies !!!
          //

          // Note: Gamma_gen should already include the trap occupancies !!!
          const double gamma_generation = viennashe::models::srh::gamma_generation(*trap_it, device, el, conf,
                                                                                   quantities, quan.get_carrier_type_id(),
                                                                                   occupancy, index_H);
          write_boundary(rhs, std::size_t(row_index),
                          - gamma_generation * trap_density * Z * volume_contribution * energy_height,
                          coupling_matrix_00,
                          viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                        );
        } // for traps

      } //assemble_traps_coupling_on_cell


      /** @brief Assembles the coupling of the distribution function with the trap occupancy for either a facet (odd-order f_lm) */
      template <typename DeviceType,
                typename TimeStepQuantitiesT,
                typename SHEQuantity,
                typename MatrixType,
                typename VectorType,
                typename CouplingMatrixType>
      void assemble_traps_coupling_on_facet(DeviceType const & device,
                                              TimeStepQuantitiesT const & quantities,
                                              SHEQuantity const & quan,
                                              viennashe::config const & conf,
                                              typename DeviceType::facet_type const & el,
                                              std::size_t index_H,
                                              MatrixType & matrix, VectorType & rhs,
                                              CouplingMatrixType const & diagonal_coupling_matrix,
                                              CouplingMatrixType const & coupling_matrix_00
                                             )
      {
        typedef typename DeviceType::mesh_type MeshType;
        typedef typename DeviceType::cell_type CellType;
        typedef typename DeviceType::trap_level_container_type       TrapContainerType;
        typedef typename TrapContainerType::const_iterator           TrapIterator;
        typedef typename viennagrid::result_of::point<MeshType>::type   PointType;
        typedef typename viennagrid::result_of::const_coboundary_range<MeshType, typename DeviceType::facet_type, CellType>::type     CellOnFacetContainer;

        (void)rhs; (void)coupling_matrix_00;
        const long row_index = quan.get_unknown_index(el, index_H);

        if (row_index < 0)
          return;

        const double Z                  = averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), el, index_H);
        const long   expansion_order    = static_cast<long>(quan.get_expansion_order(el, index_H));
        const double energy_height      = box_height(quan, el, index_H);

        viennashe::math::harmonics_iteration_type harmonics_it_id = viennashe::math::ODD_HARMONICS_ITERATION_ID;

        double volume_contribution = viennagrid::volume(el);
        volume_contribution *= detail::cell_connection_length(device.mesh(), el, viennagrid::cells(device.mesh())[0]) / static_cast<double>(PointType::dim);

        CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), el));

        for (std::size_t cell_index = 0; cell_index < cells_on_facet.size(); ++cell_index)
        {
          CellType const & cell = cells_on_facet[cell_index];

          TrapContainerType const & traps = device.get_trap_levels(cell);
          const std::size_t num_trap_unknowns    = quantities.num_trap_unknown_indices(cell);

          if (num_trap_unknowns <= 0)
            return;

          if (num_trap_unknowns != traps.size())
            throw viennashe::invalid_value_exception("The number of traps configured in the device does not match the number of unknowns for traps!", static_cast<double>(num_trap_unknowns));

          std::size_t inner_index = 0;
          for ( TrapIterator trap_it = traps.begin();
                             trap_it != traps.end();
                           ++trap_it, ++inner_index)
          {
            log::debug<log_traps>() << el << " => " << *trap_it << std::endl;

            const double occupancy = quantities.trap_occupancy(cell, inner_index);

            //
            // Term      + Gamma_rec N_T f                       (electrons)
            //             Gamma_rec N_T f                       (holes)
            //
            // Note: Gamma_rec should already include the trap occupancies !!!
            //

            const double trap_density = trap_it->density();
            // Note: Gamma_rec should already include the trap occupancies !!!
            const double gamma_recombination = viennashe::models::srh::gamma_recombination(*trap_it, device, el, conf,
                                                                                            quantities, quan.get_carrier_type_id(), occupancy, index_H);

            viennashe::util::add_block_matrix(matrix,
                                              std::size_t(row_index), std::size_t(row_index),
                                              gamma_recombination * trap_density * Z * volume_contribution * energy_height,
                                              diagonal_coupling_matrix,
                                              viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id),
                                              viennashe::math::spherical_harmonics_iterator(expansion_order, harmonics_it_id)
                                            );
          } // for traps
        }

      } // assemble_traps_coupling_on_facet

    } // namespace detail

    /** @brief Interface function for the assembly of traps
     *
     * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quantities     Container of all the SHE quantities
     * @param conf           The simulator configuration
     * @param quan           Quantity describing the trap occupancy
     * @param matrix         System matrix
     * @param rhs            Load vector
    */
    template <typename DeviceType,
              typename TimeStepQuantitiesT,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType>
    void assemble_traps_coupling(DeviceType const & device,
                                 TimeStepQuantitiesT const & quantities,
                                 viennashe::config const & conf,
                                 SHEQuantity const & quan,
                                 MatrixType & matrix,
                                 VectorType & rhs)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type        CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type       FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type          FacetIterator;

      typedef viennashe::math::sparse_matrix<double>   CouplingMatrixType;

      MeshType const & mesh = device.mesh();

      //
      // Set up scatter matrices:
      //
      std::size_t L_max = static_cast<std::size_t>(conf.max_expansion_order());
      std::size_t num_harmonics = static_cast<std::size_t>(L_max+1) * static_cast<std::size_t>(L_max+1);
      CouplingMatrixType diagonal_coupling_matrix(num_harmonics, num_harmonics);
      CouplingMatrixType coupling_matrix_00(num_harmonics, num_harmonics);
      viennashe::math::SphericalHarmonic Y_00(0, 0);

      for (std::size_t i=0; i< (L_max+1) * (L_max+1); ++i)
        diagonal_coupling_matrix(i,i) += 1.0;
      coupling_matrix_00(0,0) = 1.0 / Y_00(0,0);

        //
        // Step 1: assemble on even nodes:
        //
        CellContainer cells(mesh);
        for (CellIterator cit = cells.begin();
             cit != cells.end();
             ++cit)
        {
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            detail::assemble_traps_coupling_on_cell(device, quantities, quan, conf, *cit, index_H, matrix, rhs, diagonal_coupling_matrix, coupling_matrix_00);
        }


        //
        // Step 3: assemble on odd nodes
        //
        FacetContainer facets(mesh);
        for (FacetIterator fit = facets.begin();
             fit != facets.end();
             ++fit)
        {
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            detail::assemble_traps_coupling_on_facet(device, quantities, quan, conf, *fit, index_H, matrix, rhs, diagonal_coupling_matrix, coupling_matrix_00);
        }

    } //assemble_traps_coupling




    template <typename DeviceType,
              typename TimeStepQuantitiesT,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType>
    void assemble_traps_solver(DeviceType const & device,
                               TimeStepQuantitiesT & quantities,
                               viennashe::config const & conf,
                               SHEQuantity const & quan,
                               MatrixType & matrix,
                               VectorType & rhs
                              )
    {
      typedef typename DeviceType::mesh_type           MeshType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type        CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

      typedef typename DeviceType::trap_level_container_type     TrapContainerType;
      typedef typename TrapContainerType::const_iterator         TrapIterator;

      (void)quan; (void)matrix; (void)rhs; //eliminate unused parameter warnings

      MeshType const & mesh = device.mesh();

      CellContainer cells(mesh);
      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        const std::size_t num_trap_unknowns = quantities.num_trap_unknown_indices(*cit);
        TrapContainerType const & traps = device.get_trap_levels(*cit);
        if (num_trap_unknowns > 0)
        {

          if (num_trap_unknowns != traps.size())
            throw viennashe::invalid_value_exception("The number of traps configured in the device does not match the number of unknowns for traps!", static_cast<double>(num_trap_unknowns));

          std::size_t inner_index = 0;

          for (TrapIterator trap_it  = traps.begin();
                            trap_it != traps.end();
                          ++trap_it, ++inner_index)
          {
            // No need to assemble a matrix ...
            double new_ft = viennashe::models::srh::evaluate(*trap_it, device, *cit, conf, quantities);

            //std::cout << *cit << " => " << new_ft << std::endl;

            quantities.trap_occupancy(*cit, inner_index, new_ft);

            /* TODO: NEWTON
            matrix(unknown_index, unknown_index) = 1.0;
            if (simulation_cnt > 0)
            {
              rhs[unknown_index] = (electron_rec_rate + hole_gen_rate) / (electron_rec_rate + electron_gen_rate + hole_rec_rate + hole_gen_rate);  // aka. rate_in / (rate_in + rate_out)
              if (with_full_newton)
                rhs[unknown_index] -= current_guess[unknown_index];
            }
            else
              rhs[unknown_index] = trap_it->occupancy();

            ++unknown_index;
            */

          } //for trap_level
        } // if
      } //for cells
    }


    template <typename DeviceType,
              typename TimeStepQuantitiesT,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType>
    void assemble_traps(DeviceType const & device,
                        TimeStepQuantitiesT & quantities,
                        viennashe::config const & conf,
                        SHEQuantity const & quan,
                        MatrixType & matrix,
                        VectorType & rhs
                       )
    {

      if (conf.with_trap_selfconsistency())
        viennashe::she::assemble_traps_coupling(device, quantities, conf, quan, matrix, rhs);

      viennashe::she::assemble_traps_solver(device, quantities, conf, quan, matrix, rhs);
    }

  } //namespace she
} //namespace viennashe

#endif
