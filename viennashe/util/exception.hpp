#ifndef VIENNASHE_UTIL_EXCEPTION_HPP
#define VIENNASHE_UTIL_EXCEPTION_HPP

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

#include <string>
#include <stdexcept>

/** @file viennashe/util/exception.hpp
    @brief Defines all exceptions used/thrown in the viennashe/util/ namespace
*/

namespace viennashe
{
  namespace util
  {

    /** @brief Exception that is thrown if a string cannot be converted to the respective target type. */
    class string_conversion_exception : public std::runtime_error
    {
      public:
        string_conversion_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that a linear solver fails (typically with NaN due to zero inner product) */
    class linear_solver_exception: public std::runtime_error
    {
      public:
        linear_solver_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception thrown if for a range(a,b) the check a < b fails. */
    class invalid_range_bounds_exception : public std::runtime_error
    {
      public:
        invalid_range_bounds_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown if an invalid value is found */
    class invalid_value_exception : public std::runtime_error {
    public:
      invalid_value_exception(std::string const & str) : std::runtime_error(str) {}
    };


  } //namespace util
} //namespace viennashe

#endif
