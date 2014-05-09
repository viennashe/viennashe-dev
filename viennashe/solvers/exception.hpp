#ifndef VIENNASHE_SOLVERS_EXCEPTION_HPP
#define VIENNASHE_SOLVERS_EXCEPTION_HPP

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


#include <iostream>
#include <stdexcept>

/** @file viennashe/solvers/exception.hpp
    @brief Provides the exceptions used inside the viennashe::she namespace.
*/

namespace viennashe
{
  namespace solvers
  {

    /** @brief Exception for the case that an invalid solver is in use
      *
      * Note that only even expansion orders are allowed on vertices, and odd expansion orders on edges.
      */
    class invalid_linear_solver_exception : public std::runtime_error
    {
    public:
      invalid_linear_solver_exception() : std::runtime_error("* ViennaSHE: Invalid linear solver specified. Most likely, the linear solver is not compiled into the binary. Consider recompilation with support for the particular solver enabled.") {}
      virtual ~invalid_linear_solver_exception() throw() {}
    };


    /** @brief Exception for the case that an invalid nonlinear solver is in use
      *
      */
    class invalid_nonlinear_solver_exception : public std::runtime_error
    {
    public:
      invalid_nonlinear_solver_exception() : std::runtime_error("* ViennaSHE: Invalid nonlinear solver specified.") {}
      virtual ~invalid_nonlinear_solver_exception() throw() {}
    };

  }
}

#endif
