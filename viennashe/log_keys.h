#ifndef VIENNASHE_LOG_KEYS_H
#define VIENNASHE_LOG_KEYS_H

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

/** @file viennashe/log_keys.h
    @brief Defines the log keys used within the main viennashe:: namespace
*/

namespace viennashe
{

  /** @brief Configuration class for logging the SHE simulator */
  struct log_simulator { enum { enabled = true }; };

} // namespace viennashe


#endif

