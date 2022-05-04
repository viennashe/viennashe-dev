#ifndef VIENNASHE_SOLVERS_VIENNACL_PARALLEL_LINEAR_SOLVER_HPP
#define VIENNASHE_SOLVERS_VIENNACL_PARALLEL_LINEAR_SOLVER_HPP

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
#include "viennashe/solvers/config.hpp"
#include "src/solvers/log_keys.h"

#include "viennashe/log/log.hpp"
#include "src/solvers/log_keys.h"
#include "src/solvers/viennacl/serial_linear_solver.hpp"

// viennacl
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/io/matrix_market.hpp"

/** @file parallel_linear_solver.hpp
    @brief Provides bindings to a multithreaded block-ILU solver based on functionality in ViennaCL
*/

namespace viennashe
{
  namespace solvers
  {

    /** @brief Solves the provided system using an iterative solver with a block preconditioner on CPU
    *
    * @param system_matrix        The system matrix
    * @param rhs                  The load vector (right hand side)
    * @param config               The linear solver configuration object
    */
    template <typename CompressedMatrixType,
              typename VectorType>
    VectorType solve(CompressedMatrixType const & system_matrix,
                     VectorType const & rhs,
                     viennashe::solvers::linear_solver_config const & config,
                     viennashe::solvers::parallel_linear_solver_tag
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
      // parallel block-preconditioner
      //

      log::info<log_linear_solver>() << "* solve(): Computing block preconditioner (multi-threaded)... " << std::endl;
      //viennacl::linalg::ilut_tag precond_tag(config.ilut_entries(),
      //                                       config.ilut_drop_tolerance());
      viennacl::linalg::ilu0_tag precond_tag;

      typedef typename viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<NumericT>,
                                                           viennacl::linalg::ilu0_tag>::index_vector_type  IndexVectorType;

      IndexVectorType block_indices(config.block_preconditioner_boundaries());

      // make sure that there are block boundaries available:
      if (block_indices.size() == 0)
      {
        block_indices.resize(1);
        block_indices[0].first = 0;
        block_indices[0].second = A.size1();
      }

      viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<NumericT>,
                                          viennacl::linalg::ilu0_tag> block_preconditioner(A, precond_tag, 1);//block_indices);
      //log::debug<log_linear_solver>() << "Time: " << timer.get() << std::endl;

      //
      // Solve system:
      //
      log::info<log_linear_solver>() << "* solve(): Solving system (multi-threaded)... " << std::endl;
      viennacl::linalg::bicgstab_tag  solver_tag(config.tolerance(), config.max_iters());

      viennacl::vector<NumericT> vcl_result = viennacl::linalg::solve(A,
                                                                      b,
                                                                      solver_tag,
                                                                      block_preconditioner);

      //log::debug<log_linear_solver>() << "Time: " << timer.get() << std::endl;
      //log::debug<log_linear_solver>() << "Number of iterations (block ILUT): " << solver_tag.iters() << std::endl;

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

