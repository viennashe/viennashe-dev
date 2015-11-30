#ifndef VIENNASHE_SHE_LINEAR_SOLVER_HPP
#define VIENNASHE_SHE_LINEAR_SOLVER_HPP

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

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/she/eliminate.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/rescaling.hpp"
#include "viennashe/math/linalg_util.hpp"
#include "viennashe/solvers/forwards.h"
#include "viennashe/config.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/util/checks.hpp"


/** @file viennashe/she/linear_solver.hpp
    @brief Provides the linear solvers for SHE.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Writes an array of block indices */
    template <typename DeviceType,
              typename SHEQuantity>
    void fill_block_indices(DeviceType const & device,
                            SHEQuantity const & quan,
                            std::size_t system_size,
                            std::vector<std::pair<std::size_t, std::size_t> > & indices)
    {
      indices.resize(quan.get_value_H_size());

      std::size_t last_found = 0;

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));

      for (std::size_t index_H = 0; index_H < indices.size(); ++index_H)
      {
        indices[index_H].first = last_found;

        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          if (quan.get_unknown_index(*cit, index_H) > -1)
          {
            last_found = std::size_t(quan.get_unknown_index(*cit, index_H));
            indices[index_H].first = last_found;
            if (index_H > 0)
              indices[index_H-1].second = last_found;
            break;
          }
        }

        // make sure .second is always set:
        if (indices[index_H].second < indices[index_H].first)
          indices[index_H].second = indices[index_H].first;
      }

      indices[indices.size()-1].second = system_size;
    }



    /** @brief Public interface for solving the provided system of discretized SHE equations.
     *
     * @param device        The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan          The SHE quantities
     * @param configuration The simulator configuration
     * @param full_matrix        The system matrix containing even and odd unknowns
     * @param full_rhs           The load vector containing even and odd unknowns
     * @param reduced_unknowns   The number of reduced unknowns
    */
    template <typename DeviceType,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType>
    VectorType solve(DeviceType const & device,
                     SHEQuantity const & quan,
                     viennashe::config const & configuration,
                     MatrixType & full_matrix,
                     VectorType & full_rhs,
                     std::size_t reduced_unknowns)
    {
      typedef typename VectorType::value_type     NumericType;

      //
      // Eliminate odd unknowns
      //
      //log::debug<log_linear_solver>() << "Eliminate odd unknowns..." << std::endl;
      viennashe::math::sparse_matrix<NumericType>  compressed_matrix(reduced_unknowns, reduced_unknowns);
      VectorType compressed_rhs(reduced_unknowns);
      compressed_rhs.clear();
      compressed_rhs.resize(reduced_unknowns);
      viennashe::solvers::linear_solver_config conf(configuration.linear_solver());

      //
      // Eliminate odd unknowns
      //

      //log::info<log_linear_solver>() << "* solve(): Diagonalising odd unknowns... " << std::endl;
      viennashe::she::diagonalise_odd2odd_coupling_matrix(full_matrix, full_rhs, reduced_unknowns);

      //log::debug<log_linear_solver>() << "Full matrix: " << viennashe::util::sparse_to_string(full_matrix) << std::endl;
      //log::debug<log_linear_solver>() << "Full rhs: "    << full_rhs << std::endl;

      //log::info<log_linear_solver>() << "* solve(): Eliminating odd unknowns... " << std::endl;
      eliminate_odd_unknowns(full_matrix, full_rhs,
                            compressed_matrix, compressed_rhs);

      //log::debug<log_linear_solver>() << "Reduced matrix: " << viennashe::util::sparse_to_string(compressed_matrix) << std::endl;
      //log::debug<log_linear_solver>() << "Reduced rhs: " << compressed_rhs << std::endl;

      //
      // Check matrix for M-matrix structure (debug)
      // Makes sense for L=1 only!
      // Full system matrix is M-matrix only if no inelastic scattering is included.
      //
      //viennashe::util::m_matrix_check(compressed_matrix);

      std::vector<double> scaling_vector(compressed_matrix.size1(), 1.0);
      if (conf.scale())
      {
        setup_unknown_scaling(device, configuration, quan, scaling_vector);
        rescale_system(compressed_matrix, scaling_vector);
      }

      //
      // Normalize equation system (solver will be thankful)
      //
      viennashe::math::row_normalize_system(compressed_matrix, compressed_rhs);

      //log::debug<log_linear_solver>() << "Reduced matrix: " << viennashe::util::sparse_to_string(compressed_matrix) << std::endl;
      //log::debug<log_linear_solver>() << "Reduced rhs: " << compressed_rhs << std::endl;


      // set up preconditioner information:
      fill_block_indices(device, quan, compressed_rhs.size(), conf.block_preconditioner_boundaries());

      //
      // Solve equation system
      //
      VectorType compressed_result = viennashe::solvers::solve(compressed_matrix,
                                                               compressed_rhs,
                                                               conf);

      //
      // Scale unknowns back:
      //
      if (conf.scale())
      {
      for (std::size_t i=0; i<compressed_result.size(); ++i)
        compressed_result[i] *= scaling_vector[i];
      }

      //
      // Recover full solution of even and odd unknowns
      //
      VectorType she_result = recover_odd_unknowns(full_matrix, full_rhs, compressed_result);

      double relative_residual  = viennashe::math::norm_2(viennashe::math::subtract(viennashe::math::prod(full_matrix, she_result), full_rhs));
             relative_residual /= viennashe::math::norm_2(full_rhs);
      if (relative_residual > 1e-7)
        log::info<log_linear_solver>() << "* solve(): Relative linear solver residual of full system: " << relative_residual << std::endl;

      return she_result;
    }
  }
}

#endif
