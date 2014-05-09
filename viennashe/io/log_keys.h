#ifndef VIENNASHE_IO_LOG_KEYS_H
#define VIENNASHE_IO_LOG_KEYS_H

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

/** @file viennashe/io/log_keys.h
    @brief Defines the log keys used within the viennashe::io namespace
*/

namespace viennashe
{

  namespace io
  {
    /** @brief Tag class for logging inside the she_vtk_writer */
    struct log_she_vtk_writer { enum { enabled = true }; };
  }

} // namespace viennashe

#endif

