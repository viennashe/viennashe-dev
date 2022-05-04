#ifndef VIENNASHE_SOLVERS_VIENNACL_SERIAL_LINEAR_SOLVER_HPP
#define VIENNASHE_SOLVERS_VIENNACL_SERIAL_LINEAR_SOLVER_HPP

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
#include "viennashe/forwards.h"
#include "viennashe/util/checks.hpp"
#include "viennashe/math/linalg_util.hpp"
#include "viennashe/log/log.hpp"
#include "viennashe/solvers/config.hpp"
#include "src/solvers/log_keys.h"

// viennacl
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/io/matrix_market.hpp"


/** @file serial_linear_solver.hpp
    @brief Provides bindings to a serial ILUT-based solver shipped with ViennaCL
*/

namespace viennashe
{
  namespace solvers
  {
    namespace detail
    {
      template <typename MatrixT, typename NumericT>
      void copy(MatrixT const & assembled_matrix,
                viennacl::compressed_matrix<NumericT> & vcl_matrix)
      {
        viennacl::copy(assembled_matrix, vcl_matrix);
      }

      template <typename NumericT>
      void copy(viennashe::math::sparse_matrix<NumericT> const & assembled_matrix,
                viennacl::compressed_matrix<NumericT>          &       vcl_matrix)
      {
        std::size_t nonzeros = assembled_matrix.nnz();
        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(vcl_matrix.handle1(), assembled_matrix.size1() + 1);
        viennacl::backend::typesafe_host_array<unsigned int> col_buffer(vcl_matrix.handle2(), nonzeros);
        std::vector<NumericT> elements(nonzeros);

        std::size_t data_index = 0;

        for (std::size_t i  = 0;
                         i != assembled_matrix.size1();
                       ++i)
        {
          typedef typename viennashe::math::sparse_matrix<NumericT>::const_iterator2   AlongRowIterator;
          typedef typename viennashe::math::sparse_matrix<NumericT>::row_type          RowType;

          row_buffer.set(i, data_index);
          RowType const & row_i = assembled_matrix.row(i);

          for (AlongRowIterator col_it  = row_i.begin();
                                col_it != row_i.end();
                              ++col_it)
          {
            col_buffer.set(data_index, col_it->first);
            elements[data_index] = col_it->second;
            ++data_index;
          }
        }
        row_buffer.set(assembled_matrix.size1(), data_index);

        vcl_matrix.set(row_buffer.get(),
                       col_buffer.get(),
                       &elements[0],
                       assembled_matrix.size1(),
                       assembled_matrix.size2(),
                       nonzeros);
      }
    }


    /** @brief Solves the provided system using an iterative solver in a serial fashion (single-threaded execution)
    *
    * @param system_matrix        The system matrix
    * @param rhs                  The load vector (right hand side)
    * @param config               The linear solver configuration object
    */
    template <typename MatrixType,
              typename VectorType>
    VectorType solve(MatrixType const & system_matrix,
                     VectorType const & rhs,
                     viennashe::solvers::linear_solver_config const & config,
                     viennashe::solvers::serial_linear_solver_tag
                    )
    {
      typedef typename VectorType::value_type     NumericT;

      //
      // Step 1: Convert data to ViennaCL types:
      //
      viennacl::compressed_matrix<NumericT> A(system_matrix.size1(), system_matrix.size2());
      viennacl::vector<NumericT>            b(system_matrix.size1());

      viennacl::fast_copy(&(rhs[0]), &(rhs[0]) + rhs.size(), b.begin());
      detail::copy(system_matrix, A);


      //
      // Step 2: Setup preconditioner and run solver
      //
      log::info<log_linear_solver>() << "* solve(): Computing preconditioner (single-threaded)... " << std::endl;
      //viennacl::linalg::ilut_tag precond_tag(config.ilut_entries(),
      //                                        config.ilut_drop_tolerance());
      viennacl::linalg::ilu0_tag precond_tag;
      viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<NumericT> > preconditioner(A, precond_tag);

      log::info<log_linear_solver>() << "* solve(): Solving system (single-threaded)... " << std::endl;
      viennacl::linalg::bicgstab_tag  solver_tag(config.tolerance(), config.max_iters());

      //log::debug<log_linear_solver>() << "Compressed matrix: " << system_matrix << std::endl;
      //log::debug<log_linear_solver>() << "Compressed rhs: " << rhs << std::endl;
      viennacl::vector<NumericT> vcl_result = viennacl::linalg::solve(A,
                                                                       b,
                                                                       solver_tag,
                                                                       preconditioner);
      //log::debug<log_linear_solver>() << "Number of iterations (ILUT): " << solver_tag.iters() << std::endl;

      //
      // Step 3: Convert data back:
      //
      VectorType result(vcl_result.size());
      viennacl::fast_copy(vcl_result.begin(), vcl_result.end(), &(result[0]));

      viennashe::util::check_vector_for_valid_entries(result);

      //
      // As a check, compute residual:
      //
      log::info<log_linear_solver>() << "* solve(): residual: "
                << viennacl::linalg::norm_2(viennacl::linalg::prod(A, vcl_result) - b) / viennacl::linalg::norm_2(b)
                << " after " << solver_tag.iters() << " iterations." << std::endl;
      //log::debug<log_linear_solver>() << "SHE result (compressed): " << compressed_result << std::endl;

      return result;
    }


  } // solvers
} // viennashe

#endif

