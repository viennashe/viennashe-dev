#ifndef VIENNASHE_SOLVERS_LOG_KEYS_H
#define VIENNASHE_SOLVERS_LOG_KEYS_H

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

/** @file solvers/log_keys.h
    @brief Defines all the log keys used within the viennashe::solvers namespace
*/

namespace viennashe
{
  namespace solvers
  {

    /** @brief Logging tag for general linear solver action */
    struct log_linear_solver { enum { enabled = false, debug = false }; };

    /** @brief Logging tag for messages specific to the GPU block-ILUT preconditioner */
    struct log_gpu_block_ilut_precond { enum { enabled = false }; };

    /** @brief Logging tag for messages specific to the multithreaded (CPU-based) block-ILUT preconditioner */
    struct log_parallel_ilut_precond { enum { enabled = false }; };

  } //namespace solvers
} // namespace viennashe


#endif

