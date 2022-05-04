#ifndef VIENNASHE_SHE_SCATTERING_ALL_HPP
#define VIENNASHE_SHE_SCATTERING_ALL_HPP
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

// viennashe
#include "viennashe/she/scattering/acoustic_phonon_scattering.hpp"
#include "viennashe/she/scattering/optical_phonon_scattering.hpp"
#include "viennashe/she/scattering/impurity_scattering.hpp"
#include "viennashe/she/scattering/impact_ionization_scattering.hpp"
#include "viennashe/she/scattering/fixed_charge_scattering.hpp"
#include "viennashe/she/scattering/trapped_charge_scattering.hpp"

#include "viennashe/she/scattering/surface_scattering.hpp"

#ifdef VIENNASHE_USE_DISABLED_CODE
#include "viennashe/she/scattering/surface_acoustic_phonon_scattering.hpp"
#include "viennashe/she/scattering/surface_roughness_scattering.hpp"
#endif

/** @file viennashe/she/scattering/all.hpp
    @brief File to gather all scattering operators
*/

namespace viennashe
{
  namespace she
  {

    // Nothing in here .. this just marks the namespace one should use

  } //namespace she
} //namespace viennashe

#endif
