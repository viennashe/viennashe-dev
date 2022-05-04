#ifndef VIENNASHE_SHE_ASSEMBLE_ALL_HPP
#define VIENNASHE_SHE_ASSEMBLE_ALL_HPP

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

    template <typename DeviceType,
              typename TimeStepQuantitiesT,
              typename VertexT,
              typename EdgeT,
              typename MatrixType,
              typename VectorType>
    void assemble( DeviceType & device,
                   TimeStepQuantitiesT & old_quantities,
                   TimeStepQuantitiesT & quantities,
                   viennashe::config const & conf,
                   viennashe::she::unknown_she_quantity<VertexT, EdgeT> const & quan,
                   MatrixType & A,
                   VectorType & b,
                   bool use_timedependence, bool quan_valid)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

      typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
      typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

      typedef viennashe::math::sparse_matrix<double>   CouplingMatrixType;

      typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_quantity_type      SpatialUnknownType;

      std::vector< scattering_base<DeviceType> * > scattering_processes;

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
        MeshType const & mesh = device.mesh();

        SpatialUnknownType const &     potential =     quantities.get_unknown_quantity(viennashe::quantity::potential());
        SpatialUnknownType const & old_potential = old_quantities.get_unknown_quantity(viennashe::quantity::potential());  //TODO: Take old timestep

        viennashe::she::unknown_she_quantity<VertexT, EdgeT> const & old_quan = old_quantities.she_quantity(quan.get_name());

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
          scattering_processes.push_back(new acoustic_phonon_scattering<DeviceType>(device, conf));
        }

        if (conf.scattering().optical_phonon().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Optical phonon scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new optical_phonon_scattering<DeviceType>(device, conf, conf.energy_spacing()));
        }

        if (conf.scattering().ionized_impurity().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Ionized impurity scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new ionized_impurity_scattering<DeviceType>(device, conf));
        }

        if (conf.scattering().impact_ionization().enabled())
        {
          // Warn the user if we already know that he/she is going to simulate bullshit.
          if ( ! conf.with_holes() || conf.get_hole_equation() != viennashe::EQUATION_SHE )
            log::warn() << std::endl << "WARNING: II scattering enabled, but 'BTE for holes' is disabled! Expect inconsistent results!" << std::endl;
          if ( ! conf.with_electrons() || conf.get_electron_equation() != viennashe::EQUATION_SHE )
            log::warn() << std::endl << "WARNING: II scattering enabled, but 'BTE for electrons' is disabled! Expect inconsistent results!" << std::endl;

          scattering_processes.push_back(new impact_ionization_scattering<DeviceType>(device, conf));
        }

        if (conf.with_traps() && conf.scattering().trapped_charge().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Trapped charge scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new trapped_charge_scattering<DeviceType, TimeStepQuantitiesT>(device, conf, quantities));
        }

        typedef typename viennashe::electric_field_wrapper<DeviceType, SpatialUnknownType> ElectricFieldAccessor;
        ElectricFieldAccessor Efield(device, potential);

        if (conf.scattering().surface().enabled())
        {
          log::debug<log_assemble_all>() << "assemble(): Surface roughness scattering is ENABLED!" << std::endl;
          scattering_processes.push_back(new surface_scattering<DeviceType, ElectricFieldAccessor>(device, conf, Efield));
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

        CellContainer cells(mesh);
        for (CellIterator cit = cells.begin();
            cit != cells.end();
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
                                                 *cit, index_H,
                                                 scatter_op_in, scatter_op_out);
          }

          //
          // Free streaming operator L{f}
          //

          // iterate over neighbor cells holding the odd unknowns:
          FacetOnCellContainer facets_on_cell(*cit);
          for (FacetOnCellIterator focit = facets_on_cell.begin();
              focit != facets_on_cell.end();
              ++focit)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling coupling with facet " << *focit << std::endl;

            CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *focit, *cit);
            if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

            CouplingMatrixType coupling_matrix_diffusion = coupling_matrix_in_direction(a_x, a_y, a_z,
                                                                                        *cit,
                                                                                        *other_cell_ptr,
                                                                                        quan.get_carrier_type_id());

            // note that the sign change due to MEDS is included in the choice of the normal vector direction (order of vertex vs. other_vertex):
            // - B \cdot n   for even unknowns,
            // + B \cdot n   for odd unknowns
            CouplingMatrixType coupling_matrix_drift = coupling_matrix_in_direction(b_x, b_y, b_z, *other_cell_ptr, *cit, quan.get_carrier_type_id());

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
                                                      old_potential);
            }
          }
        } //for cells


        //
        // Step 2: Assemble on odd 'nodes' (i.e. facets). TODO: Resolve code duplication w.r.t. above
        //
        log::info<log_assemble_all>() << "* assemble_all(): Odd unknowns..." << std::endl;

        FacetContainer facets(mesh);
        for (FacetIterator fit = facets.begin();
             fit != facets.end();
             ++fit)
        {
          if (log_assemble_all::enabled)
            log::debug<log_assemble_all>() << "* assemble_all(): Assembling on facet " << *fit << std::endl;

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
                                                *fit, index_H,
                                                scatter_op_in, scatter_op_out);
          }

          //
          // Free streaming operator L{f}
          //

          // iterate over cells of facet
          CellOnFacetContainer cells_on_facet(mesh, fit.handle());
          for (CellOnFacetIterator cofit  = cells_on_facet.begin();
                                   cofit != cells_on_facet.end();
                                 ++cofit)
          {
            if (log_assemble_all::enabled)
              log::debug<log_assemble_all>() << "* assemble_all(): Assembling coupling with cell " << *cofit << std::endl;

            CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *fit, *cofit);
            if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

            CouplingMatrixType coupling_matrix_diffusion = coupling_matrix_in_direction(a_x_transposed, a_y_transposed, a_z_transposed,
                                                                                        *other_cell_ptr,
                                                                                        *cofit,
                                                                                        quan.get_carrier_type_id());

            // note that the sign change due to MEDS is included in the choice of the normal vector direction (order of vertex vs. other_vertex):
            // - B \cdot n   for even unknowns,
            // + B \cdot n   for odd unknowns
            CouplingMatrixType coupling_matrix_drift = coupling_matrix_in_direction(b_x_transposed, b_y_transposed, b_z_transposed, *other_cell_ptr, *cofit, quan.get_carrier_type_id());

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
                                                      old_potential);
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
