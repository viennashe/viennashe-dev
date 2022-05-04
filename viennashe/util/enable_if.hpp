#ifndef VIENNASHE_UTIL_ENABLE_IF_HPP
#define VIENNASHE_UTIL_ENABLE_IF_HPP

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

/** @file viennashe/util/enable_if.hpp
    @brief Simple enable-if variant that uses the SFINAE pattern
*/

namespace viennashe
{
    /** @brief Simple enable-if variant that uses the SFINAE pattern */
    template <bool b, class T = void>
    struct enable_if
    {
      typedef T   type;
    };

    template <class T>
    struct enable_if<false, T> {};

} //namespace viennashe


#endif
