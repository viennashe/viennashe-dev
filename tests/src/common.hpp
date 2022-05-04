#ifndef VIENNASHE_TESTS_COMMON_HPP
#define	VIENNASHE_TESTS_COMMON_HPP
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

/** @file tests/src/common.hpp
    @brief Contains common functions, functors and other classes often needed by the tests
 */

#include "viennashe/forwards.h"
#include "viennashe/log/log.hpp"

#include <cmath>
#include <algorithm>

namespace viennashe
{
  namespace testing
  {

    /**
     * @brief Performs a fuzzy (up to a tolerance) equal. Returns true if is and should are equal within the tolerance tol.
     * @param is The value
     * @param should The value is should be
     * @param tol The tolerance
     * @return True if is == should within the tolerance tol, else false.
     */
    inline bool fuzzy_equal(double is, double should, double tol = 1e-1)
    {
      if ( (is < should || is > should)
          && std::fabs(is - should) / std::max(std::abs(is), std::fabs(should)) > tol)
      {
        viennashe::log::error() << "fuzzy_equal(): tol    = " << tol    << std::endl;
        viennashe::log::error() << "fuzzy_equal(): is     = " << is     << std::endl;
        viennashe::log::error() << "fuzzy_equal(): should = " << should << std::endl;
        viennashe::log::error() << "fuzzy_equal(): diff   = " << std::fabs(is - should) / std::max(std::fabs(is), std::fabs(should)) << std::endl;
        return false;
      }
      return true;
    }

  } // namespace tests
} // namespace viennashe

#endif	/* VIENNASHE_TESTS_COMMON_HPP */

