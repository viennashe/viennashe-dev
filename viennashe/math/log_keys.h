#ifndef VIENNASHE_MATH_LOG_KEYS_H
#define VIENNASHE_MATH_LOG_KEYS_H

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

/** @file viennashe/math/log_keys.h
    @brief Defines the log keys used within the viennashe::math namespace
*/

namespace viennashe
{
  namespace math
  {
    /** @brief Controls logging for the derivative of spherical harmonics with respect to theta */
    struct log_SphericalHarmonic_dTheta { enum { enabled = false }; };

    /** @brief Controls logging within the associated Legendre polynomial class */
    struct log_AssocLegendre { enum { enabled = false }; };

    /** @brief Controls logging within the numerical integration routine */
    struct log_integrate { enum { enabled = false }; };
  }
}

#endif

