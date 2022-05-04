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

#include <iostream>
#include <cstdlib>
#include <vector>

#ifndef NDEBUG
#define NDEBUG
#endif

#include "viennashe/forwards.h"

#include "src/solvers/solve.hpp"
#include "viennashe/solvers/forwards.h"
#include "viennashe/math/linalg_util.hpp"

#include "src/solvers/native/dense_linear_solver.hpp"
#include "src/solvers/petsc/petsc_solver.hpp"

namespace viennashe
{
  namespace solvers
  {

    // sparse matrix
    std::vector<double> solve(viennashe::math::sparse_matrix<double> & A,
        std::vector<double> const & b, linear_solver_config const & config)
    {
      return solve_impl(A, b, config);
    }

    // dense matrix
    std::vector<double> solve(viennashe::math::dense_matrix<double> & A,
        std::vector<double> const & b)
    {
      linear_solver_config dummy_config;
      return viennashe::solvers::solve(A, b, dummy_config,
          viennashe::solvers::dense_linear_solver_tag());
    }

//    // PETSC
//    std::vector<double>
//       solve(viennashe::math::sparse_matrix<double> & A,
//             std::vector<double> const & b,
//             linear_solver_config const & config)
//       {
//         linear_solver_config dummy_config(linear_solver_ids::petsc_parallel_linear_solver); //
//         return  solve_impl(A, b, dummy_config);
//       }

  }

}
