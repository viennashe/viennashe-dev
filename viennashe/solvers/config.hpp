#ifndef VIENNASHE_SOLVERS_CONFIG_HPP
#define VIENNASHE_SOLVERS_CONFIG_HPP

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

#include <iostream>
#include <algorithm>

#include "viennashe/forwards.h"

#include "viennashe/solvers/exception.hpp"

/** @file viennashe/solvers/config.hpp
 @brief The SHE configuration class is defined here.
 */

namespace viennashe
{
  namespace solvers
  {

    //
    // Part 1: Linear solvers
    //

    /** @brief Provides IDs for the linear solvers */
    struct linear_solver_ids
    {
        enum
        {
          dense_linear_solver, /// Gauss solver (use for systems up to ~1000-by-1000 only)
          serial_linear_solver,        /// single-threaded solver (using ILUT)
          parallel_linear_solver,      /// multi-threaded solver (block ILUT)
          gpu_parallel_linear_solver, /// gpu-assisted solver (block ILUT)
          petsc_parallel_linear_solver, /// PETSC-assisted solver (Hypre AMG)
	  petsc_parallel_AMGX_solver   /// gpu/PETSC-assisted solver
        };
    };

    /** @brief A configuration class holding options for use within the different linear solvers
     *
     *  Internal note: May require a bit more polishing, particularly with respect to preconditioners
     */
    class linear_solver_config: public linear_solver_ids
    {
      public:
        typedef std::vector<std::pair<std::size_t, std::size_t> > block_preconditioner_boundaries_container;

        linear_solver_config()
            : id_(linear_solver_ids::serial_linear_solver), tol_(1e-13), max_iters_(
                1000), ilut_entries_(60), ilut_drop_tol_(1e-4), do_scale_(true)
        {
        }

        /** @brief Controls the scaling of unkowns. Give "false" to disable unkown scaling */
        void scale(bool value)
        {
          do_scale_ = value;
        }
        /** @brief Returns whether or not the unkowns are scaled before solving the equation set */
        bool scale() const
        {
          return do_scale_;
        }

        /** @brief Returns the current linear solver ID */
        long id() const
        {
          return id_;
        }
        /** @brief Sets a new linear solver. Provide IDs defined in @see linear_solver_ids */
        void set(long solver_id)
        {
          switch (solver_id)
          {
            case dense_linear_solver:
            case serial_linear_solver:
            case parallel_linear_solver:
            case petsc_parallel_linear_solver:
            case petsc_parallel_AMGX_solver:
#ifdef VIENNASHE_HAVE_GPU_SOLVER
              case gpu_parallel_linear_solver:
#endif
              id_ = solver_id;
              break;
            default:
              throw invalid_linear_solver_exception();
          }
        }

        /** @brief Sets a new linear solver. Provide IDs defined in @see linear_solver_ids */
        void set(std::string solver_name)
        {
          // to upper to ease the comparison and to be case insensitive
          std::transform(solver_name.begin(), solver_name.end(),
              solver_name.begin(), ::tolower);

          if (solver_name == "dense_linear_solver")
          {
            id_ = dense_linear_solver;
          } else if (solver_name == "serial_linear_solver")
          {
            id_ = serial_linear_solver;
          } else if (solver_name == "parallel_linear_solver")
          {
            id_ = parallel_linear_solver;
          } else if (solver_name == "petsc_parallel_linear_solver")
          {
            id_ = petsc_parallel_linear_solver;
          }
          else if (solver_name == "petsc_parallel_AMGX_solver")
          {
            id_ = petsc_parallel_AMGX_solver;
          }
#ifdef VIENNASHE_HAVE_GPU_SOLVER
          else if (solver_name == "gpu_parallel_linear_solver")
          {
            id_ = gpu_parallel_linear_solver;
          }
#endif
          else
          {
            throw invalid_linear_solver_exception();
          }
        }
        void finalize_petsc() {


        }
        /** @brief Returns the tolerance for the (iterative) linear solver */
        double tolerance() const
        {
          return tol_;
        }
        /** @brief Sets the tolerance for the (iterative) linear solver. Non-positive values are ignored */
        void tolerance(double tol)
        {
          if (tol > 0)
            tol_ = tol;
        }

        /** @brief Returns the maximum number of iterative linear solver iterations */
        std::size_t max_iters() const
        {
          return max_iters_;
        }
        /** @brief Sets the maximum number of iterations within the iterative linear solver */
        void max_iters(std::size_t iters)
        {
          max_iters_ = iters;
        }

        /** @brief Returns the maximum number of entries in each of the two factors of the ILUT preconditioner */
        std::size_t ilut_entries() const
        {
          return ilut_entries_;
        }
        /** @brief Sets the maximum number of entries in each of the two factors of the ILUT preconditioner */
        void ilut_entries(std::size_t num)
        {
          ilut_entries_ = num;
        }

        /** @brief Returns the drop tolerance for the ILUT preconditioenr */
        double ilut_drop_tolerance() const
        {
          return ilut_drop_tol_;
        }
        /** @brief Sets the drop tolerance for the ILUT preconditioner. Non-positive values are ignored */
        void ilut_drop_tolerance(double tol)
        {
          if (tol > 0)
            ilut_drop_tol_ = tol;
        }

        block_preconditioner_boundaries_container & block_preconditioner_boundaries()
        {
          return block_precond_boundaries_;
        }
        block_preconditioner_boundaries_container const & block_preconditioner_boundaries() const
        {
          return block_precond_boundaries_;
        }
        
        int getArgc() const
        {
          return argc;
        }
        
        void setArgc(int argc)
        {
          this->argc = argc;
        }
        
        char** getArgv() const
        {
          return argv;
        }
        
        void setArgv(char** argv)
        {
          this->argv = argv;
        }

        void finalize() const
        {
        //  viennashe::solvers::PETSC_handler<NumericT,VectorType>->finalize();
        }
      private:
        //linear solver:
        long id_;
        double tol_;
        int argc;
        char **argv;
        std::size_t max_iters_;
        std::size_t ilut_entries_;
        double ilut_drop_tol_;
        block_preconditioner_boundaries_container block_precond_boundaries_;
        bool do_scale_;
    };

    //
    // Part 2: Nonlinear solvers
    //

    /** @brief Provides IDs for the nonlinear iteration schemes */
    struct nonlinear_solver_ids
    {
        enum
        {
          gummel_nonlinear_solver,      /// Gummel iteration
          newton_nonlinear_solver       /// Newton's method
        };
    };

    /** @brief A configuration class holding options for use within the different linear solvers
     *
     *  Internal note: May require a bit more polishing, particularly with respect to preconditioners
     */
    class nonlinear_solver_config: public nonlinear_solver_ids
    {
      public:
        nonlinear_solver_config()
            : id_(nonlinear_solver_ids::gummel_nonlinear_solver), iterN(0), THRESHOLD(400), tol_(1e-8), max_iters_(
                100), damping_(0.3)
        {
        }

        /** @brief Returns the current linear solver ID */
        long id() const
        {
          return
              iterN > THRESHOLD ?
                  nonlinear_solver_ids::newton_nonlinear_solver : id_;
          //return id_;
        }
        /** @brief Sets a new linear solver. Provide IDs defined in @see linear_solver_ids */
        void set(long solver_id)
        {
          switch (solver_id)
          {
            case nonlinear_solver_ids::gummel_nonlinear_solver:
            case nonlinear_solver_ids::newton_nonlinear_solver:
              id_ = solver_id;
              break;
            default:
              throw invalid_nonlinear_solver_exception();
          }
        }

        /** @brief Sets a new nonlinear solver. */
        void set(std::string solver_name)
        {
          // to upper to ease the comparison and to be case insensitive
          std::transform(solver_name.begin(), solver_name.end(),
              solver_name.begin(), ::tolower);

          if (solver_name == "gummel_nonlinear_solver"
              || solver_name == "gummel")
          {
            id_ = nonlinear_solver_ids::gummel_nonlinear_solver;
          } else if (solver_name == "newton_nonlinear_solver"
              || solver_name == "newton")
          {
            id_ = nonlinear_solver_ids::newton_nonlinear_solver;
          } else
          {
            throw invalid_nonlinear_solver_exception();
          }
        }

        /** @brief Returns the tolerance for the nonlinear solver potential update */
        double tolerance() const
        {
          return tol_;
        }
        /** @brief Sets the tolerance for the nonlinear solver potential update. Non-positive values are ignored. */
        void tolerance(double tol)
        {
          if (tol > 0)
            tol_ = tol;
        }

        /** @brief Returns the maximum number of nonlinear solver iterations */
        std::size_t max_iters() const
        {
          return max_iters_;
        }
        /** @brief Sets the maximum number of iterations within the nonlinear solver */
        void max_iters(std::size_t iters)
        {
          max_iters_ = iters;
        }

        /** @brief Returns the damping for the nonlinear solver */
        double damping() const
        {
          return damping_;
        }
        /** @brief Sets the tolerance for the nonlinear solver. Negative values are ignored. */
        void damping(double d)
        {
          if (damping_ >= 0)
            damping_ = d;
        }

        /** @brief Increment the Number of Iterations for the nonlinear solver.  */

        void operator++(int)
        {
          iterN++;
        }

        long getThreshold() const
        {
          return THRESHOLD;
        }

        void setThreshold(long threshold)
        {
          THRESHOLD = threshold;
        }

        ;

      private:
        long id_;
        long iterN;
        long THRESHOLD;
        double tol_;
        std::size_t max_iters_;
        double damping_;
    };

  }
}

#endif
