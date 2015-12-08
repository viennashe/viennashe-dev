#ifndef VIENNASHE_SHE_ASSEMBLE_ALL_HPP
#define VIENNASHE_SHE_ASSEMBLE_ALL_HPP

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

#include "viennashe/she/assemble_common.hpp"
#include "viennashe/she/assemble_boundary.hpp"
#include "viennashe/she/assemble_streaming.hpp"
#include "viennashe/she/assemble_scattering.hpp"
#include "viennashe/she/assemble_timederivative.hpp"
#include "viennashe/she/scattering/assemble_ee_scattering.hpp"
#include "viennashe/she/assemble_traps.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"

#include "viennashe/postproc/electric_field.hpp"

/** @file viennashe/she/assemble_all.hpp
    @brief Generic assembly of traps is implemented here
*/

namespace viennashe
{
  namespace she
  {

    template <typename TimeStepQuantitiesT,
              typename MatrixType,
              typename VectorType>
    void assemble( viennashe::device & device,
                   TimeStepQuantitiesT & old_quantities,
                   TimeStepQuantitiesT & quantities,
                   viennashe::config const & conf,
                   viennashe::she::unknown_she_quantity<double> const & quan,
                   MatrixType & A,
                   VectorType & b,
                   bool use_timedependence, bool quan_valid)
    {
      typedef viennashe::math::sparse_matrix<double>   CouplingMatrixType;

      typedef typename TimeStepQuantitiesT::unknown_quantity_type      SpatialUnknownType;

      std::vector< scattering_base * > scattering_processes;

      if (conf.with_traps())
      {
        if (! conf.with_electrons() || ! conf.with_holes())
          throw viennashe::unavailable_feature_exception("Trapping without considering electrons or holes is not supported!");
        if ( conf.get_electron_equation() != viennashe::EQUATION_SHE)
          throw viennashe::unavailable_feature_exception("Trapping without SHE for electrons is not supported!");
        if ( conf.get_hole_equation() != viennashe::EQUATION_SHE)
          throw viennashe::unavailable_feature_exception("Trapping without SHE for holes is not supported!");
      }

//      try
//      {
        viennagrid_mesh mesh = device.mesh();

        SpatialUnknownType const &     potential =     quantities.get_unknown_quantity(viennashe::quantity::potential());
        SpatialUnknownType const & old_potential = old_quantities.get_unknown_quantity(viennashe::quantity::potential());  //TODO: Take old timestep

        viennashe::she::unknown_she_quantity<double> const & old_quan = old_quantities.she_quantity(quan.get_name());

        //
        // Set up scatter matrices:
        //
        const std::size_t L_max = static_cast<std::size_t>(conf.max_expansion_order());
        const std::size_t num_harmonics = std::size_t(L_max+1) * std::size_t(L_max+1);
        CouplingMatrixType scatter_op_in(num_harmonics, num_harmonics);
        CouplingMatrixType scatter_op_out(num_harmonics, num_harmonics);

        for (std::size_t i=0; i < std::size_t(L_max+1) * std::size_t(L_max+1); ++i)
          scatter_op_out(i,i) += 1.0;
        scatter_op_in(0,0) += 1.0;

        //// preprocessing: compute coefficients a_{l,m}^{l',m'} and b_{l,m}^{l',m'}
        std::size_t Lmax = static_cast<std::size_t>(conf.max_expansion_order());            //maximum expansion order
        std::size_t coupling_rows = static_cast<std::size_t>((Lmax+1) * (Lmax+1));
        std::size_t coupling_cols = coupling_rows;

        log::debug<log_assemble_all>() << "* assemble_all(): Computing coupling matrices..." << std::endl;
        CouplingMatrixType identity(coupling_rows, coupling_cols);
        for (std::size_t i=0; i<coupling_rows; ++i)
          for (std::size_t j=0; j<coupling_cols; ++j)
            identity(i,j) = (i == j) ? 1.0 : 0.0;

        CouplingMatrixType a_x(coupling_rows, coupling_cols);
        CouplingMatrixType a_y(coupling_rows, coupling_cols);
        CouplingMatrixType a_z(coupling_rows, coupling_cols);


        CouplingMatrixType b_x(coupling_rows, coupling_cols);
        CouplingMatrixType b_y(coupling_rows, coupling_cols);
        CouplingMatrixType b_z(coupling_rows, coupling_cols);

        //note: interchanged coordinates
        fill_coupling_matrices(a_x, a_y, a_z,
                               b_x, b_y, b_z,
                               static_cast<int>(Lmax));

        CouplingMatrixType a_x_transposed = a_x.trans();
        CouplingMatrixType a_y_transposed = a_y.trans();
        CouplingMatrixType a_z_transposed = a_z.trans();

        CouplingMatrixType b_x_transposed = b_x.trans();
        CouplingMatrixType b_y_transposed = b_y.trans();
        CouplingMatrixType b_z_transposed = b_z.trans();

        if (log_assemble_all::enabled && log_assemble_all::debug)
        {
          log::debug<log_assemble_all>() << "a_x: " << a_x << std::endl;
          log::debug<log_assemble_all>() << "a_y: " << a_y << std::endl;
          log::debug<log_assemble_all>() << "a_z: " << a_z << std::endl;
          log::debug<log_assemble_all>() << "b_x: " << b_x << std::endl;
          log::debug<log_assemble_all>() << "b_y: " << b_y << std::endl;
          log::debug<log_assemble_all>() << "b_z: " << b_z << std::endl;

          log::debug<log_assemble_all>() << "identity: " << identity << std::endl;

          log::debug<log_assemble_all>() << "scatter_op_out: " << scatter_op_out << std::endl;
          log::debug<log_assemble_all>() << "scatter_op_in: "  << scatter_op_in << std::endl;
        }

        //
        // Setup vector of scattering processes:
        //


        if (conf.scattering().acoustic_phonon().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Acoustic phonon scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new acoustic_phonon_scattering(device, conf));
        }

        if (conf.scattering().optical_phonon().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Optical phonon scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new optical_phonon_scattering(device, conf, conf.energy_spacing()));
        }

        if (conf.scattering().ionized_impurity().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Ionized impurity scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new ionized_impurity_scattering(device, conf));
        }

        if (conf.scattering().impact_ionization().enabled())
        {
          // Warn the user if we already know that he/she is going to simulate bullshit.
          if ( ! conf.with_holes() || conf.get_hole_equation() != viennashe::EQUATION_SHE )
            log::warn() << std::endl << "WARNING: II scattering enabled, but 'BTE for holes' is disabled! Expect inconsistent results!" << std::endl;
          if ( ! conf.with_electrons() || conf.get_electron_equation() != viennashe::EQUATION_SHE )
            log::warn() << std::endl << "WARNING: II scattering enabled, but 'BTE for electrons' is disabled! Expect inconsistent results!" << std::endl;

          scattering_processes.push_back(new impact_ionization_scattering(device, conf));
        }

        if (conf.with_traps() && conf.scattering().trapped_charge().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Trapped charge scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new trapped_charge_scattering<TimeStepQuantitiesT>(device, conf, quantities));
        }

        typedef typename viennashe::electric_field_wrapper<SpatialUnknownType> ElectricFieldAccessor;
        ElectricFieldAccessor Efield(device, potential);

        if (conf.scattering().surface().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Surface roughness scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new surface_scattering<ElectricFieldAccessor>(device, conf, Efield));
        }


        //
        // Assemble SHE system:
        //    - scattering operators on vertices
        //    - free streaming operator on vertices
        //    - scattering operators on edges
        //    - free streaming operator on edges
        //    - any other stuff (traps on cells, etc.)
        //

        if (quan_valid && conf.scattering().electron_electron() && conf.with_electrons())
        {
          log::debug<log_assemble_all>() << "assemble(): Electron electron scattering is ENABLED!" << std::endl;
          assemble_ee_scattering(device, conf, quan, old_quan, A, b);
        }

        //
        // Step 1: Assemble on even nodes
        //
        log::debug<log_assemble_all>() << "* assemble_all(): Even unknowns..." << std::endl;

        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

        viennagrid_element_id *cells_begin, *cells_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          log::debug<log_assemble_all>() << "* assemble_all(): Assembling on cell " << *cit << std::endl;

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling on energy index " << index_H << std::endl;

            assemble_boundary_on_box(device, conf, quan,
                                     A, b,
                                     *cit, index_H,
                                     identity);

            if (viennashe::materials::is_conductor(device.get_material(*cit)))
              continue;

            //
            // Scattering operator Q{f}
            //
            assemble_scattering_operator_on_box( scattering_processes,
                                                 device, conf, quan,
                                                 A, b,
                                                 *cit, index_H, 1.0,
                                                 scatter_op_in, scatter_op_out,
                                                 false);
          }

          //
          // Free streaming operator L{f}
          //

          // iterate over neighbor cells holding the odd unknowns:
          viennagrid_element_id *facets_on_cell_begin, *facets_on_cell_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(mesh, *cit, cell_dim - 1, &facets_on_cell_begin, &facets_on_cell_end));
          for (viennagrid_element_id *focit  = facets_on_cell_begin;
                                      focit != facets_on_cell_end;
                                    ++focit)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling coupling with facet " << *focit << std::endl;

            viennagrid_element_id *other_cell_ptr;
            util::get_other_cell_of_facet(mesh, *focit, *cit, &other_cell_ptr);
            if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

            CouplingMatrixType coupling_matrix_diffusion = coupling_matrix_in_direction(a_x, a_y, a_z,
                                                                                        mesh, *cit, *other_cell_ptr,
                                                                                        quan.get_carrier_type_id());

            // note that the sign change due to MEDS is included in the choice of the normal vector direction (order of vertex vs. other_vertex):
            // - B \cdot n   for even unknowns,
            // + B \cdot n   for odd unknowns
            CouplingMatrixType coupling_matrix_drift = coupling_matrix_in_direction(b_x, b_y, b_z,
                                                                                    mesh, *other_cell_ptr, *cit,
                                                                                    quan.get_carrier_type_id());

            for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            {
              assemble_free_streaming_operator_on_box( device, conf, quan,
                                                       A, b,
                                                       *cit, *focit, index_H,
                                                       coupling_matrix_diffusion,
                                                       coupling_matrix_drift,
                                                       false);
            }
          } //for edges

          //
          // Time dependence df/dt (and possibly df/dH * dH/dt)
          //
          if (use_timedependence)
          {
            for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            {
              viennashe::she::assemble_timederivative(device, conf, quan, old_quan,
                                                      A, b,
                                                      *cit, index_H,
                                                      identity,
                                                      potential,
                                                      old_potential, false);
            }
          }
        } //for cells


        //
        // Step 2: Assemble on odd 'nodes' (i.e. facets). TODO: Resolve code duplication w.r.t. above
        //
        log::info<log_assemble_all>() << "* assemble_all(): Odd unknowns..." << std::endl;

        viennagrid_element_id *facets_begin, *facets_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim-1, &facets_begin, &facets_end));

        for (viennagrid_element_id *fit  = facets_begin;
                                    fit != facets_end;
                                  ++fit)
        {
          if (log_assemble_all::enabled)
            log::debug<log_assemble_all>() << "* assemble_all(): Assembling on facet " << *fit << std::endl;

          viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

          if (cells_on_facet_begin + 1 == cells_on_facet_end) // facet is on boundary, hence no contribution
            continue;

          viennagrid_numeric cell_connection_len;
          std::vector<viennagrid_numeric> centroid_1(3);
          std::vector<viennagrid_numeric> centroid_2(3);
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, cells_on_facet_begin[0], &(centroid_1[0])));
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, cells_on_facet_begin[1], &(centroid_2[0])));
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(cell_dim, &(centroid_1[0]), &(centroid_2[0]), &cell_connection_len));

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling on energy index " << index_H << std::endl;

            //
            // Scattering operator Q{f}
            //
            assemble_scattering_operator_on_box(scattering_processes,
                                                device, conf, quan,
                                                A, b,
                                                *fit, index_H, cell_connection_len / double(cell_dim),
                                                scatter_op_in, scatter_op_out, true);
          }

          //
          // Free streaming operator L{f}
          //

          // iterate over cells of facet
          for (viennagrid_element_id *cofit  = cells_on_facet_begin;
                                      cofit != cells_on_facet_end;
                                    ++cofit)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling coupling with cell " << *cofit << std::endl;

            viennagrid_element_id other_cell = (cells_on_facet_begin[0] == *cofit) ? cells_on_facet_begin[1] : cells_on_facet_begin[0];

            CouplingMatrixType coupling_matrix_diffusion = coupling_matrix_in_direction(a_x_transposed, a_y_transposed, a_z_transposed,
                                                                                        mesh, other_cell, *cofit,
                                                                                        quan.get_carrier_type_id());

            // note that the sign change due to MEDS is included in the choice of the normal vector direction (order of vertex vs. other_vertex):
            // - B \cdot n   for even unknowns,
            // + B \cdot n   for odd unknowns
            CouplingMatrixType coupling_matrix_drift = coupling_matrix_in_direction(b_x_transposed, b_y_transposed, b_z_transposed,
                                                                                    mesh, other_cell, *cofit,
                                                                                    quan.get_carrier_type_id());

            for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            {
              assemble_free_streaming_operator_on_box( device, conf, quan,
                                                       A, b,
                                                       *cofit, *fit, index_H,
                                                       coupling_matrix_diffusion,
                                                       coupling_matrix_drift,
                                                       true);
            }
          } //for vertices

          //
          // Time dependence df/dt (and possibly df/dH * dH/dt)
          //
          if (use_timedependence)
          {
            for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
            {
              viennashe::she::assemble_timederivative(device, conf, quan, old_quan,
                                                      A, b,
                                                      *fit, index_H,
                                                      identity,
                                                      potential,
                                                      old_potential,
                                                      true);
            }
          }

        } //for facets


        // Assemble traps on cells (to be integrated into the assembly above):
        if (conf.with_traps())
        {
          log::debug<log_assemble_all>() << "assemble(): Assembly for traps ..." << std::endl;
          viennashe::she::assemble_traps(device, quantities, conf, quan, A, b);
        }

        //
        // Cleanup:
        //
        for (std::size_t i=0; i<scattering_processes.size(); ++i)
        {
          if ( scattering_processes[i] ) delete scattering_processes[i];
          scattering_processes[i] = 0;
        }
/*      }
      catch (...)
      {
        //
        // Cleanup:
        //
        for (std::size_t i=0; i<scattering_processes.size(); ++i)
          if ( scattering_processes[i] ) delete scattering_processes[i];
        // Rethrow
        throw;
      }
*/

    }
  } //namespace she
} //namespace viennashe

#endif
