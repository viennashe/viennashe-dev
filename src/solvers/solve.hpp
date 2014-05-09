#ifndef VIENNASHE_SOLVERS_SOLVE_HPP
#define VIENNASHE_SOLVERS_SOLVE_HPP

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

// std
#include <stdexcept>

#include "viennashe/solvers/config.hpp"
#include "viennashe/solvers/exception.hpp"
#include "src/solvers/viennacl/all.h"
#include "src/solvers/native/dense_linear_solver.hpp"

#include "viennashe/util/checks.hpp"
#include "viennacl/linalg/norm_2.hpp"

/** @file solve.hpp
    @brief Provides a generic solver layer, providing bindings to a dense Gauss solver, ViennaCL, etc.
*/

namespace viennashe
{
  /** @brief Namespace containing a variety of different linear solvers. */
  namespace solvers
  {


    /** @brief Public interface for solving a system of linear equations
    *
    * @param system_matrix     The system matrix containing even and odd unknowns
    * @param rhs               The load vector containing even and odd unknowns
    * @param config            Linear solver configuration object
    */
    template <typename MatrixType,
              typename VectorType>
    VectorType solve_impl(MatrixType & system_matrix,
                          VectorType const & rhs,
                          linear_solver_config const & config)
    {
      // check for invalid entries first
      const long invalid_row = viennashe::util::matrix_consistency_check(system_matrix);
      if (invalid_row >= 0)
      {
        throw viennashe::util::linear_solver_exception("solve(): Found empty row in system_matrix.");
      }
      viennashe::util::check_vector_for_valid_entries(rhs, "solve(): ");
      //check for trivial solution (some solvers have problems with that, hence we do it explicitly here)
      if (!viennacl::linalg::norm_2(rhs))
      {
        VectorType result(rhs.size());
        std::fill(result.begin(), result.end(), 0);
        return result;
      }

      switch (config.id())
      {
        case linear_solver_config::dense_linear_solver:
          return viennashe::solvers::solve(system_matrix, rhs, config, viennashe::solvers::dense_linear_solver_tag());
        case linear_solver_config::serial_linear_solver:
          return viennashe::solvers::solve(system_matrix, rhs, config, viennashe::solvers::serial_linear_solver_tag());
        case linear_solver_config::parallel_linear_solver:
          return viennashe::solvers::solve(system_matrix, rhs, config, viennashe::solvers::parallel_linear_solver_tag());
#ifdef VIENNASHE_HAVE_GPU_SOLVER
        case linear_solver_config::gpu_parallel_linear_solver:
          return solve(system_matrix, rhs, config, viennashe::solvers::gpu_parallel_linear_solver_tag());
#endif
        default:
          throw invalid_linear_solver_exception();
      }
    }

  }
}

#endif
