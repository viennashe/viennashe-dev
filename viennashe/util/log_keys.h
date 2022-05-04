#ifndef VIENNASHE_UTIL_LOG_KEYS_H
#define VIENNASHE_UTIL_LOG_KEYS_H

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

/** @file viennashe/util/log_keys.h
    @brief Defines the log keys used within the viennashe::util namespace
*/

namespace viennashe
{
  namespace util
  {
    /** @brief Logging key for messages occurring during M-matrix check */
    struct log_m_matrix_check { enum { enabled = true }; };

    /** @brief Logging key for messages from the built-in ortho-grid generator */
    struct log_generate_device { enum { enabled = true }; };
  }
}

#endif

