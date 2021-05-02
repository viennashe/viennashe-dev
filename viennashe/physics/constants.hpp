#ifndef VIENNASHE_PHYSICS_CONSTANTS_HPP
#define VIENNASHE_PHYSICS_CONSTANTS_HPP

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

//
// ALL UNITS ARE SI !
//


/** @file viennashe/physics/constants.hpp
    @brief Provides a number of fundamental constants. All constants in SI units.
*/

namespace viennashe
{
  /** @brief Namespace for physics models used within ViennaSHE */
  namespace physics
  {
    /** @brief Implementations. Not intended to be used by a library user. */
    namespace detail
    {
      /** @brief Implementation class holding basic physics constants */
      template <bool dummy = true>  //template argument in order to control linkage
      struct constants
      {
        /** @brief Speed of light in vacuum */
        static const double c_0;
        /** @brief Electron rest mass */
        static const double mass_electron;
        /** @brief Elementary charge */
        static const double q;
        /** @brief Boltzmann constant */
        static const double kB;
        /** @brief Planck constant */
        static const double h;
        /** @brief Modified Planck constant */
        static const double hbar;
        /** @brief Permittivity of vacuum. */
        static const double eps_0;
        /** @brief Permeability of vacuum. */
        static const double mu_0;
        /** @brief Intrinsic carrier concentration in silicon. */
        static const double ni;
      };

      template <bool b>
      const double constants<b>::c_0 = 299792458.0; // m / s
      template <bool b>
      const double constants<b>::mass_electron = 9.10938291e-31; // kg
      template <bool b>
      const double constants<b>::q      = 1.602176565e-19;     // As
      template <bool b>
      const double constants<b>::kB     = 1.3806488e-23;    // J/K
      template <bool b>
      const double constants<b>::h      = 6.62606896e-34; // Js
      template <bool b>
      const double constants<b>::hbar   = 6.62606896e-34 / (2.0 * 3.1415926535897932384626433832795); // Js
      template <bool b>
      const double constants<b>::eps_0  =  1.0e7 / ( 4.0 * 3.1415926535897932384626433832795 * 299792458.0 * 299792458.0) ; // = 10^7 / (4*pi*c_0^2) ~ 8.854e-12 As / Vm
      template <bool b>
      const double constants<b>::mu_0  =  ( 4.0 * 3.1415926535897932384626433832795 ) * 1.0e-7 ; // = 10^-7 * (4*pi) ~ 12.56e-7   Vs / Am
      template <bool b>
      const double constants<b>::ni = 1.08e16; // m^-3
    }

    typedef detail::constants<>   constants;

  } //namespace physics
} //namespace viennashe
#endif
