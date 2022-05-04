#ifndef SRC_SOLVERS_NATIVE_PETSC_SOLVER_HPP_
#define SRC_SOLVERS_NATIVE_PETSC_SOLVER_HPP_

/* ============================================================================
 Copyright (c) 2011-2019, Institute for Microelectronics,
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
#include "viennashe/util/timer.hpp"
#include <petscksp.h>
#include <mpi.h>
#include "petsc.hpp"

namespace viennashe
  {
    namespace solvers
      {
        template <typename NumericT, typename VectorType>
          VectorType
          solve (viennashe::math::sparse_matrix<NumericT> const &A,
              VectorType const &b,
              viennashe::solvers::linear_solver_config const &config,
              viennashe::solvers::petsc_linear_solver_tag)
          {

            viennashe::util::timer stopwatch;

            typedef viennashe::solvers::PETSC_matrix<NumericT, VectorType> PETSC_MATRIX;
            typedef viennashe::solvers::PETSC_vector<NumericT, VectorType> PETSC_VECTOR;
            typedef viennashe::solvers::PETSC_handler<NumericT, VectorType> PETSC_HANDLER;
            PETSC_HANDLER &ref = PETSC_HANDLER::get_handler (config.getArgc (),
                                                             config.getArgv ());
            // Check size nnz


            // copy A and b to work system Ax = b
            PETSC_MATRIX Ap (A); // work matrix
            PetscBarrier (NULL);

            //printf ("\n Matric Pointer: %p \n",&Ap );
            PETSC_VECTOR Bp (b, Ap);
            PetscBarrier (NULL);

            stopwatch.start ();
            Ap.solve (Bp);
            PetscBarrier (NULL);
            KSPConvergedReason reason = Ap.get_Creason ();
            if (reason < 0)
              {
                log::info<log_linear_solver> () << "Divergence.\n";
              }
            else
              {
                log::info<log_linear_solver> () << "Convergence in"
                    << Ap.get_Citer () << "iterations." << std::endl;
              }

            VectorType newVector(Bp);
            printf ("\n [%d]Solver Time (%d): %f Reason %d\n",ref.getRankW(),Ap.get_Citer() ,stopwatch.get (),reason);
            return newVector;
          }

      }
  }
#endif /* SRC_SOLVERS_NATIVE_PETSC_SOLVER_HPP_ */
