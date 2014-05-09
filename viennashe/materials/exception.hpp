#ifndef VIENNASHE_MATERIALS_EXCEPTION_HPP
#define VIENNASHE_MATERIALS_EXCEPTION_HPP

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

// std
#include <stdexcept>
#include <string>


/** @file viennashe/materials/exception.hpp
    @brief Contains all the exceptions used for accessing materials.
*/

namespace viennashe
{
  namespace materials
  {
    /** @brief Exception for the case that an invalid material is in use */
    class invalid_material_exception: public std::runtime_error
    {
      public:
        invalid_material_exception(std::string const & str) : std::runtime_error(str) {}
    };

  } //namespace materials
} //namespace viennashe

#endif

