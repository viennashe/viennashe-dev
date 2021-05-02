#ifndef VIENNASHE_MATH_CONSTANTS_HPP
#define VIENNASHE_MATH_CONSTANTS_HPP

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

/** @file viennashe/math/constants.hpp
    @brief Provides a number of fundamental math constants.
*/

namespace viennashe
{
  namespace math
  {
    /** @brief Implementations. Not intended to be used by a library user. */
    namespace detail
    {
      /** @brief Implementation class holding basic math constants */
      template <bool dummy = true>  //template argument in order to control linkage
      struct constants
      {
        /** @brief Pi */
        static const double pi;

        /** @brief Euler's number */
        static const double e;

      };

      template <bool b>
      const double constants<b>::pi = 3.1415926535897932384626433832795; // 1

      template <bool b>
      const double constants<b>::e  = 2.718281828459045235360287471352; // 1
    }

    /** @brief Convenience typedef for accessing mathematical constants, e.g. viennashe::math::constants::pi */
    typedef detail::constants<>   constants;

  } //namespace math
} //namespace viennashe
#endif
