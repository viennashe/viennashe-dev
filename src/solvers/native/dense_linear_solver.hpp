#ifndef VIENNASHE_SOLVERS_NATIVE_DENSE_LINEAR_SOLVER_HPP
#define VIENNASHE_SOLVERS_NATIVE_DENSE_LINEAR_SOLVER_HPP

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

// viennashe
#include "viennashe/util/checks.hpp"
#include "viennashe/math/linalg_util.hpp"
#include "viennashe/log/log.hpp"
#include "viennashe/solvers/config.hpp"
#include "src/solvers/log_keys.h"


/** @file dense_linear_solver.hpp
    @brief Implements a dense Gauss solver
*/


namespace viennashe
{
  namespace solvers
  {

    namespace detail
    {
      // pivoting strategy: Pick element of B(:,j) below diagonal with largest modulus
      template <typename NumericT, typename VectorType>
      void pivot_matrix(viennashe::math::dense_matrix<NumericT> & B,
                        VectorType & c,
                        std::vector<std::size_t> pivot_vector,
                        std::size_t current_row)
      {
        // find pivot row:
        std::size_t pivot_row = current_row;
        NumericT    pivot_element = B(current_row, current_row);

        for (std::size_t i=current_row; i<B.size1(); ++i)
        {
          double element = B(i, current_row);
          if (std::fabs(element) > std::fabs(pivot_element))
          {
            pivot_row = i;
            pivot_element = element;
          }
        }

        // swap rows:
        if (pivot_row > current_row)
        {
          std::swap(pivot_vector[pivot_row], pivot_vector[current_row]);
          std::swap(c[pivot_row], c[current_row]);
          for (std::size_t j=0; j<B.size2(); ++j)
            std::swap(B(pivot_row, j), B(current_row, j));
        }
      }
    }

    /** @brief Solves the provided system using a dense Gauss solver with partial pivoting.
    *
    * @param A        The system matrix
    * @param b        The load vector (right hand side)
    * @param config   The linear solver configuration object
    */
    template <typename NumericT, typename VectorType>
    VectorType solve(viennashe::math::dense_matrix<NumericT> const & A,
                     VectorType const & b,
                     viennashe::solvers::linear_solver_config const & config,
                     viennashe::solvers::dense_linear_solver_tag)
    {
      (void)config; //Silence unused parameter warnings
      log::info<log_linear_solver>() << "* solve(): Solving system (Gauss solver, single-threaded)... " << std::endl;

      // copy A and b to work system Bx = c
      viennashe::math::dense_matrix<NumericT> B(A); // work matrix
      VectorType c(b);
      std::vector<std::size_t> pivot_vector(B.size1());

      for (std::size_t i=0; i<pivot_vector.size(); ++i)
        pivot_vector[i] = i;

      //
      // Phase 1: Eliminate lower-triangular entries to obtain upper triangular matrix:
      //
      for (std::size_t i=0; i<B.size1(); ++i)
      {
        // pivot if necessary:
        detail::pivot_matrix(B, c, pivot_vector, i);

        if (!B(i,i))
          throw std::runtime_error("Provided matrix is singular!");

        // eliminate column entries below diagonal
        for (std::size_t elim_row = i+1; elim_row < B.size1(); ++elim_row)
        {
          NumericT factor = B(elim_row, i) / B(i,i);
          c[elim_row] -= c[i] * factor;
          for (std::size_t j=i; j<B.size2(); ++j)
            B(elim_row, j) -= B(i, j) * factor;
        }
      }

      //
      // Phase 2: Substitute:
      //
      for (std::size_t i=0; i < B.size1(); ++i)
      {
        std::size_t row = B.size1() - (i+1);
        for (std::size_t j=row+1; j<B.size2(); ++j)
          c[row] -= B(row, j) * c[j];
        c[row] /= B(row, row);
      }

      //
      // Phase 3: Form actual solution vector by inverting the pivoting:
      //
      VectorType result(c.size());
      for (std::size_t i=0; i<result.size(); ++i)
        result[pivot_vector[i]] = c[i];

      return result;
    }


    /** @brief Specialization for a sparse matrix */
    template <typename NumericT,
              typename VectorType>
    VectorType solve(viennashe::math::sparse_matrix<NumericT> const & A,
                     VectorType const & b,
                     viennashe::solvers::linear_solver_config const & config,
                     viennashe::solvers::dense_linear_solver_tag)
    {
      // convert to dense matrix and solve
      viennashe::math::dense_matrix<NumericT> A_dense(A.size1(), A.size2());

      for (std::size_t i  = 0;
                       i != A.size1();
                     ++i)
      {
        typedef typename viennashe::math::sparse_matrix<NumericT>::const_iterator2   AlongRowIterator;
        typedef typename viennashe::math::sparse_matrix<NumericT>::row_type          RowType;

        RowType const & row_i = A.row(i);

        for (AlongRowIterator col_it  = row_i.begin();
                              col_it != row_i.end();
                            ++col_it)
        {
          A_dense(i, col_it->first) = col_it->second;
        }
      }

      return solve(A_dense, b, config, viennashe::solvers::dense_linear_solver_tag());
    }

  } // namespace solvers
} // namespace viennashe


#endif
