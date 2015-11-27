#ifndef VIENNASHE_SHE_ASSEMBLE_STREAMING_HPP
#define VIENNASHE_SHE_ASSEMBLE_STREAMING_HPP

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
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/config.hpp"
#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/assemble_common.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/util/filter.hpp"

#include "viennashe/log/log.hpp"

/** @file viennashe/she/assemble_streaming.hpp
    @brief Assembly of the free-streaming operator is implemented here
*/

namespace viennashe
{
  namespace she
  {

    inline double compute_weighted_interface_area(viennagrid_mesh mesh,
                                                  viennagrid_element_id cell,
                                                  viennagrid_element_id other_cell,
                                                  viennagrid_element_id facet)
    {
      std::vector<double> centroid_1(3), centroid_2(3);

      double connection_len;
      viennagrid_element_centroid(mesh, cell,       &(centroid_1[0]));
      viennagrid_element_centroid(mesh, other_cell, &(centroid_2[0]));

      for (std::size_t i=0; i<centroid_1.size(); ++i)
        centroid_1[i] -= centroid_2[i];
      viennagrid_norm_2(3, &(centroid_1[0]), &connection_len);

      std::vector<double> cell_connection_normalized(3);
      for (std::size_t i=0; i<cell_connection_normalized.size(); ++i)
        cell_connection_normalized[i] = centroid_1[i] / connection_len;

      std::vector<double> facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(mesh, cell, facet);

      double facet_volume;
      viennagrid_element_volume(mesh, facet, &facet_volume);
      double normal_share;
      viennagrid_inner_prod(3, &(facet_unit_normal[0]), &(cell_connection_normalized[0]), &normal_share);

      return facet_volume * normal_share;
    }

    /** @brief Worker function for the assembly of the free streaming operator. Handles the assembly of both even and odd unknowns.
     *
     * @tparam HarmonicsIter1    Iterator over the indices of the spherical harmonics to be assembled. Even harmonics for even unknowns, odd harmonics for odd unknowns
     * @tparam HarmonicsIter2    Iterator over the indices of the coupled spherical harmonics. Opposite parity of HarmonicsIter1 due to even-to-odd and odd-to-even property of the free streaming operator
     * @tparam DeviceType        The device descriptor class
     * @tparam SHEController     The SHE controller type giving access to all SHE-specific quantities (unknown indices, etc.)
     * @tparam MatrixType        The system matrix type
     * @tparam VectorType        The load vector type
     * @tparam ElementType1      The topological element type the assembly is associated with. Vertex for even unknowns, edge for odd unknowns.
     * @tparam ElementType2      The topological element type of coupled unknowns. Edge for even unknowns, vertex for odd unknowns. Remember even-to-odd (odd-to-even) property of free streaming operator.
     * @tparam CouplingMatrixType Type of the coupling matrices a_{l,m}^{l',m'} and b_{l,m}^{l',m'}
     */
    template <typename DeviceType,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType,
              typename CouplingMatrixType>
    void assemble_free_streaming_operator_on_box(DeviceType const & device,
                                                 viennashe::config const & conf,
                                                 SHEQuantity const & quan,
                                                 MatrixType & A,
                                                 VectorType & /*b*/,
                                                 viennagrid_element_id cell, viennagrid_element_id facet, std::size_t index_H,
                                                 CouplingMatrixType const & coupling_matrix_diffusion,
                                                 CouplingMatrixType const & coupling_matrix_drift,
                                                 bool odd_assembly)
    {
      typename viennashe::config::dispersion_relation_type dispersion = conf.dispersion_relation(quan.get_carrier_type_id());

      //const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

      viennagrid_element_id *other_cell_ptr;
      util::get_other_cell_of_facet(device.mesh(), facet, cell, &other_cell_ptr);

      if (!other_cell_ptr) return;

      long row_index = odd_assembly ? quan.get_unknown_index(facet, index_H) : quan.get_unknown_index(cell, index_H);

      if (row_index < 0) //no unknown here
        return;

      long col_index =  odd_assembly ? quan.get_unknown_index(cell, index_H) : quan.get_unknown_index(facet, index_H);

      if (col_index < 0) //other element does not carry an unknown, so nothing to do here
        return;

      const double weighted_interface_area = compute_weighted_interface_area(device.mesh(), cell, *other_cell_ptr, facet);

      long expansion_order_row    = static_cast<long>(odd_assembly ? quan.get_expansion_order(facet, index_H) : quan.get_expansion_order(cell,  index_H));
      long expansion_order_column = static_cast<long>(odd_assembly ? quan.get_expansion_order(cell,  index_H) : quan.get_expansion_order(facet, index_H));

      double lower_kinetic_energy = lower_kinetic_energy_at_facet(device, quan, facet, index_H);
      double upper_kinetic_energy = upper_kinetic_energy_at_facet(device, quan, facet, index_H);

      viennashe::math::harmonics_iteration_type harmonics_it_id_row = (odd_assembly ? viennashe::math::ODD_HARMONICS_ITERATION_ID
                                                                                    : viennashe::math::EVEN_HARMONICS_ITERATION_ID);
      viennashe::math::harmonics_iteration_type harmonics_it_id_col = (odd_assembly ? viennashe::math::EVEN_HARMONICS_ITERATION_ID
                                                                                    : viennashe::math::ODD_HARMONICS_ITERATION_ID);

      //
      // 'diffusion' term
      //

      double matrix_coeff_diffusion = 0;
      switch (conf.she_discretization_type())
      {
      case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
        matrix_coeff_diffusion = integral_vZ(quan, dispersion, facet, index_H) * weighted_interface_area;
        break;
      case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
        matrix_coeff_diffusion = integral_v(quan, dispersion, facet, index_H) * weighted_interface_area;
        break;
      default: throw std::runtime_error("assemble_free_streaming_operator_on_box(): Unknown SHE discretization type!");
      }

      // write matrix contribution:
      if (col_index < 0)
      {
        log::error() << "* assemble_free_streaming_operator(): Invalid column index detected while assembling diffusion operator. " << std::endl;
        log::error<log_assemble_free_streaming_operator>() << " Diagnostics: " << std::endl << "index_H: " << index_H << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "cell: " << cell << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "row_index: " << row_index << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "facet: " << facet << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "odd_assembly: " << odd_assembly << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "DOF el2: " << col_index << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "other_cell: " << *other_cell_ptr << std::endl;
      }

      viennashe::util::add_block_matrix(A, std::size_t(row_index), std::size_t(col_index),
                                        matrix_coeff_diffusion,
                                        coupling_matrix_diffusion,
                                        viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                        viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col)
                                      );

      if (log_assemble_free_streaming_operator::debug)
      {
        log::debug<log_assemble_free_streaming_operator>() << "* assemble_free_streaming_operator_on_box(): diffusion term contribution:" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "      int_vZ * interface_area" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "  at (" << row_index << ", " << col_index << ")" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << " matrix_coeff_diffusion = " << matrix_coeff_diffusion << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "weighted_interface_area = " << weighted_interface_area << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << " coupling_matrix_diffusion: " << coupling_matrix_diffusion << std::endl;
      }

      /*
      if (with_full_newton)
      {
        double potential_offset = 0.05;
        double energy_offset = potential_offset * viennashe::physics::constants::q;
        double dint_vZ_dphi = (integral_vZ(lower_kinetic_energy + energy_offset,
                                           upper_kinetic_energy + energy_offset,
                                           dispersion)
                                - integral_vZ(lower_kinetic_energy,
                                              upper_kinetic_energy,
                                              dispersion)
                              ) / potential_offset;

        double d_dphi = dint_vZ_dphi * interface_area;

        if (odd_assembly)
        {
          // coupling with 'considered vertex' (edge-to-vertex coupling) is negative for electrons (energy, negative for holes
          if (potential_index(vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                     row_index, potential_index(vertex),
                                                     (assemble_electrons ? -0.5 : 0.5) * d_dphi,
                                                     coupling_matrix_diffusion,
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col),
                                                     current_guess, col_index);

          // coupling with 'other' vertex (other vertex of edge) is negative for electrons, positive for holes:
          if (potential_index(other_vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                     row_index, potential_index(other_vertex),
                                                     (assemble_electrons ? 0.5 : -0.5) * d_dphi,
                                                     coupling_matrix_diffusion,
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col),
                                                     current_guess, col_index);
        }
        else
        {
          // coupling with 'this' vertex is negative for electrons, positive for holes
          if (potential_index(vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                     row_index, potential_index(vertex),
                                                     (assemble_electrons ? -0.5 : 0.5) * d_dphi,
                                                     coupling_matrix_diffusion,
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col),
                                                     current_guess, col_index);

          // coupling with 'other' vertex is positive for electrons, negative for holes:
          if (potential_index(other_vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                     row_index, potential_index(other_vertex),
                                                     (assemble_electrons ? 0.5 : -0.5) * d_dphi,
                                                     coupling_matrix_diffusion,
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                                     viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col),
                                                     current_guess, col_index);
        }

        // residual contribution:
        viennashe::util::subtract_folded_block_vector(b,
                                                      row_index,
                                                      matrix_coeff_diffusion,
                                                      coupling_matrix_diffusion,
                                                      viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                                      viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col),
                                                      current_guess, col_index);
      }
      */

      //
      // 'drift' term
      //

      double connection_len;
      {
        std::vector<double> centroid_1(3), centroid_2(3);
        viennagrid_element_centroid(device.mesh(), cell,            &(centroid_1[0]));
        viennagrid_element_centroid(device.mesh(), *other_cell_ptr, &(centroid_2[0]));

        for (std::size_t i=0; i<centroid_1.size(); ++i)
          centroid_1[i] -= centroid_2[i];
        viennagrid_norm_2(3, &(centroid_1[0]), &connection_len);
      }


      double Z             = averaged_density_of_states(quan, dispersion, facet, index_H);
      double int_Z_over_hk = integral_Z_over_hk(lower_kinetic_energy,
                                                upper_kinetic_energy,
                                                dispersion);
      double force = force_on_facet_accessor()(quan, device.mesh(), *other_cell_ptr, cell, index_H); //TODO: Here we need the global force, not just the projection on the edge
      //force_along_edge(controller, other_vertex, vertex, index_H);

      double matrix_coeff_drift = 0;
      switch (conf.she_discretization_type())
      {
      case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
        matrix_coeff_drift = int_Z_over_hk * force * weighted_interface_area * connection_len / 2.0;
        break;
      case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
        matrix_coeff_drift = int_Z_over_hk / Z * force * weighted_interface_area * connection_len / 2.0;
        break;
      default: throw std::runtime_error("assemble_free_streaming_operator_on_box(): Unknown SHE discretization type!");
      }


      if (col_index < 0)
      {
        log::error() << "* assemble_free_streaming_operator(): Invalid column index detected while assembling drift operator. Diagnostics: " << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "index_H: " << index_H << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "cell: " << cell << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "row_index: " << row_index << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "facet: " << facet << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "odd_assembly: " << odd_assembly << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "DOF el2: " << col_index << std::endl;
        log::error<log_assemble_free_streaming_operator>() << "other_cell: " << *other_cell_ptr << std::endl;
      }

      viennashe::util::add_block_matrix(A, std::size_t(row_index), std::size_t(col_index),
                                        matrix_coeff_drift,
                                        coupling_matrix_drift,
                                        viennashe::math::spherical_harmonics_iterator(expansion_order_row,    harmonics_it_id_row),
                                        viennashe::math::spherical_harmonics_iterator(expansion_order_column, harmonics_it_id_col)
                                       );

      if (log_assemble_free_streaming_operator::debug)
      {
        log::debug<log_assemble_free_streaming_operator>() << "* assemble_free_streaming_operator_on_box(): drift term contribution:" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "      int_Z_over_hk / Z * force * interface_area * edge_len / 2.0" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "  at (" << row_index << ", " << col_index << ")" << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "    int_Z_over_hk = " << int_Z_over_hk << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "            force = " << force << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << " weighted_interface_area = " << weighted_interface_area << std::endl;
        log::debug<log_assemble_free_streaming_operator>() << "   connection_len = " << connection_len << std::endl;
      }


      /*if (with_full_newton)
      {
        //We assume 'rate' to be independent of the potential, which is true for phonon scattering

        double potential_offset = 0.05;
        double energy_offset = potential_offset * viennashe::physics::constants::q;

        double dint_Z_over_hk_dphi = (integral_Z_over_hk(lower_kinetic_energy + energy_offset, upper_kinetic_energy + energy_offset, dispersion)
                                      - integral_Z_over_hk(lower_kinetic_energy, upper_kinetic_energy, dispersion)
                                     ) / potential_offset;

        double d_dphi =   dint_Z_over_hk_dphi * force * interface_area * edge_len / 2.0; //Note: Ignoring dependence of force on the potential for the moment (term does not contribute for first-order SHE anyways)

        std::vector<double> fold_vector(coupling_matrix_drift.size2());
        if (assemble_electrons)
          current_guess.fill(viennashe::ELECTRON_TYPE_ID, el2, controller.kinetic_energy(el2, index_H, viennashe::ELECTRON_TYPE_ID), index_H, fold_vector);
        else
          current_guess.fill(viennashe::HOLE_TYPE_ID,     el2, controller.kinetic_energy(el2, index_H, viennashe::HOLE_TYPE_ID),     index_H, fold_vector);

        if (odd_assembly)
        {
          // coupling with 'considered vertex' (edge-to-vertex coupling) is positive for electrons, negative for holes
          if (potential_index(vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                    row_index, potential_index(vertex),
                                                    (assemble_electrons ? 0.5 : -0.5) * d_dphi,
                                                    coupling_matrix_drift,
                                                    HarmonicsIter1(expansion_order_row),
                                                    HarmonicsIter2(expansion_order_column),
                                                    fold_vector);

          // coupling with 'other' vertex (other vertex of edge) is negative for electrons, positive for holes:
          if (potential_index(other_vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                    row_index, potential_index(other_vertex),
                                                    (assemble_electrons ? -0.5 : 0.5) * d_dphi,
                                                    coupling_matrix_drift,
                                                    HarmonicsIter1(expansion_order_row),
                                                    HarmonicsIter2(expansion_order_column),
                                                    fold_vector);
        }
        else
        {
          // coupling with 'this' vertex is negative for electrons, positive for holes
          if (potential_index(vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                    row_index, potential_index(vertex),
                                                    (assemble_electrons ? -1.0 : 1.0) * d_dphi,
                                                    coupling_matrix_drift,
                                                    HarmonicsIter1(expansion_order_row),
                                                    HarmonicsIter2(expansion_order_column),
                                                    fold_vector);

          // coupling with 'other' vertex is positive for electrons, negative for holes:
          if (potential_index(other_vertex) >= 0)
            viennashe::util::add_folded_block_matrix(A,
                                                    row_index, potential_index(other_vertex),
                                                    (assemble_electrons ? 1.0 : -1.0) * d_dphi,
                                                    coupling_matrix_drift,
                                                    HarmonicsIter1(expansion_order_row),
                                                    HarmonicsIter2(expansion_order_column),
                                                    fold_vector);
        }

        // residual contribution:
        viennashe::util::subtract_folded_block_vector(b,
                                                      row_index,
                                                      matrix_coeff,
                                                      coupling_matrix_drift,
                                                      HarmonicsIter1(expansion_order_row),
                                                      HarmonicsIter2(expansion_order_column),
                                                      fold_vector);
      }*/

    }



  } //namespace she
} //namespace viennashe

#endif
