#ifndef VIENNASHE_SOLVERS_FORWARDS_H
#define VIENNASHE_SOLVERS_FORWARDS_H

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

#include <vector>
#include "viennashe/math/linalg_util.hpp"

#include "viennashe/solvers/config.hpp"

/** @file viennashe/solvers/forwards.h
    @brief Forward declarations for a generic solver layer, providing bindings to a native Gauss solver, ViennaCL, etc.
*/

namespace viennashe
{
  /** @brief Namespace containing a variety of different linear solvers. */
  namespace solvers
  {
    /** @brief Public interface for solving a system of linear equations represented using a sparse matrix
    *
    * @param A        The system matrix containing even and odd unknowns
    * @param b        The load vector containing even and odd unknowns
    * @param config   Linear solver configuration object
    */
    std::vector<double>
    solve(viennashe::math::sparse_matrix<double> & A,
          std::vector<double> const & b,
          linear_solver_config const & config);


    /** @brief Public interface for solving a system of linear equations using a dense matrix
    *
    * @param A        The system matrix containing even and odd unknowns
    * @param b        The load vector containing even and odd unknowns
    */
    std::vector<double>
    solve(viennashe::math::dense_matrix<double> & A,
          std::vector<double> const & b);
  }
}

#endif
