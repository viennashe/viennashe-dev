#ifndef VIENNASHE_MATH_EXCEPTION_HPP
#define VIENNASHE_MATH_EXCEPTION_HPP

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
#include <iostream>
#include <fstream>
#include <stdexcept>

/** @file viennashe/math/exception.hpp
    @brief All the exceptions used within the viennashe::math namespace.
*/

namespace viennashe
{
  /** @brief Namespace for all math specific code  */
  namespace math
  {

    /** @brief Exception which is thrown if the size of a solution vector is truncated (too small or big) */
    class truncated_solution_vector_size_exception : public std::runtime_error {
    public:
      truncated_solution_vector_size_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception which is thrown if an iterate cannot be retrieved from the vector of iterates (unknowns) */
    class iterate_not_found_exception : public std::runtime_error {
    public:
      iterate_not_found_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown if a (linear) solver provided to some nonlinear solver is not valid (e.g. Nullpointer provided) */
    class provided_solver_invalid_exception : public std::runtime_error {
    public:
      provided_solver_invalid_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown if an updater for a nonlinear solver is not valid (e.g. Nullpointer provided) */
    class provided_updater_invalid_exception : public std::runtime_error {
    public:
      provided_updater_invalid_exception(std::string const & str) : std::runtime_error(str) {}
    };

  } //namespace math
}//namespace viennashe


#endif

