#ifndef VIENNASHE_MODELS_ALL_HPP
#define VIENNASHE_MODELS_ALL_HPP

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

// viennashe
#include "viennashe/forwards.h"

// exceptions
#include "viennashe/models/exception.hpp"

// models
#include "viennashe/models/bandgap.hpp"
#include "viennashe/models/carrier_mass.hpp"
#include "viennashe/models/mobility.hpp"
#include "viennashe/models/hot_carrier.hpp"

//#include "viennashe/models/random_dopant_model.hpp"

// ...
// ... add your model here ...
// ...


/** @file viennashe/models/all.hpp
    @brief Contains the basic models interface.
 */

namespace viennashe
{
  /** @brief Namespace containing everything models related. */
  namespace models
  {

    // MODEL REGISTRY HERE

  } // models
} // viennashe

#endif /* VIENNASHE_MODELS_ALL_HPP */

