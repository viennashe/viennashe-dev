#ifndef VIENNASHE_PHYSICS_PHYSICS_HPP
#define VIENNASHE_PHYSICS_PHYSICS_HPP

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
#include <math.h>
#include <cmath>

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/math/constants.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/exception.hpp"
#include "viennashe/materials/all.hpp"


/** @file viennashe/physics/physics.hpp
    @brief Returns a few helper routines for computing physical quantities. To be replaced in the future.
 */

namespace viennashe
{
  namespace physics
  {

    /** @brief Auxiliary namespace for converting quantities of different units
     *
     */
    namespace convert
    {

      //
      // Replace this by a sane unit framework!
      //

      inline double joule_to_eV(const double eps)
      {
        return eps / viennashe::physics::constants::q;
      }

      inline double eV_to_joule(const double eps)
      {
        return eps * viennashe::physics::constants::q;
      }
    }

    /** @brief Returns the thermal potential for the provided temperature */
    inline double get_thermal_potential(double T)
    {
      return (constants::kB * T) / constants::q;
    }

    /** @brief Returns the optical phonon number for use with optical phonon scattering
     *
     * @param eps_optical    The phonon energy
     * @param temperature    The lattice temperature
     */
    inline double get_optical_phonon_number(double eps_optical, double temperature = 300.0 /* Kelvin */)
    {
      const double kB = constants::kB;
      return 1.0 / (std::exp(eps_optical / (kB * temperature)) - 1.0);
    }

    /**
     * @brief Returns the band weight (N_c or N_v)
     * @param temperature The lattice temperature
     * @param ctype Determines whether to return N_c (electrons) or N_v (holes)
     * @return The band weight
     *
     */
    inline double get_band_weight(double temperature, viennashe::carrier_type_id ctype)
    {
      const double Ncv0 = 2.50933670314359713002e19 * 100 * 100 * 100;
      // Insert carrier mass model here
      const double me = (ctype == viennashe::ELECTRON_TYPE_ID) ? 0.33 : 0.6;
      return Ncv0 * std::pow(me * (temperature / 300.0), 1.5) * 6.0;
    }

    /**
     * @brief Returns the band edge relative to the reference energy (mid-gap)
     * @param ctype Give electrons for E_c and holes for E_v
     * @return The band edge in Joule
     *
     */
    inline double get_band_edge(viennashe::carrier_type_id const & ctype)
    {
      return ((ctype == viennashe::ELECTRON_TYPE_ID) ? 1.0 : -1.0) * viennashe::materials::si::band_gap() * 0.5;
    }


    /**
     * @brief Returns the thermal velocity at the given lattice temperature for a given carrier type
     * @param temperature The lattice temperature
     * @param ctype The charge carrier type
     * @return The thermal velocity in meters per second
     */
    inline double get_thermal_velocity(double temperature, viennashe::carrier_type_id ctype)
    {
      // TODO: use carrier mass model
      (void)ctype; // fix unused parameter warning
      const double kBT = viennashe::physics::constants::kB * temperature;
      return std::sqrt(8.0 * kBT / (viennashe::math::constants::pi * viennashe::physics::constants::mass_electron) );
    }

    /**
     * @brief Returns the auxilary carrier concentration
     * @param temperature The lattice temperature
     * @param ctype The carrier type (electrons or holes)
     * @return A carrier concentration
     */
    inline double get_auxilary_concentration(double temperature, viennashe::carrier_type_id ctype)
    {
      const double kbT  = viennashe::physics::constants::kB * temperature;
      const double band_edge   = viennashe::physics::get_band_edge(ctype);
      const double band_weight = viennashe::physics::get_band_weight(temperature, ctype);
      const double polarity = (ctype == viennashe::ELECTRON_TYPE_ID) ? -1.0 : 1.0;

      return band_weight * std::exp(polarity * band_edge / kbT);
    }

    /**
     * @brief Computes the built-in potential for a given temperature and given doping
     * @param temperature The lattice temperature
     * @param doping_n The donor doping (SI units only!)
     * @param doping_p The acceptor doping (SI units only!)
     * @return The built-in potential in Volt
     */
    inline double built_in_potential(double temperature, double doping_n, double doping_p)
    {
      if (doping_n <= 0 && doping_p <= 0) return 0;

      const double net_doping = (doping_n - doping_p);
      const double T    = temperature;
      const double naux = viennashe::physics::get_auxilary_concentration(T, viennashe::ELECTRON_TYPE_ID);
      const double paux = viennashe::physics::get_auxilary_concentration(T, viennashe::HOLE_TYPE_ID);

      double x = std::sqrt(net_doping * net_doping + 4.0 * naux * paux);

      if ( net_doping >= 0.0 )
      {
        x = (x + net_doping) * 0.5 / naux;
        return  T * std::log(x) * viennashe::physics::constants::kB / viennashe::physics::constants::q;
      }
      else
      {
        x = (x - net_doping) * 0.5 / paux;
        return -T * std::log(x) * viennashe::physics::constants::kB / viennashe::physics::constants::q;
      }

    }

    /**
     * @brief Consistently calculates the contact carrier concentrations for thermal bath contacts
     * @param temperature The lattice temperature
     * @param doping_n The donor doping
     * @param doping_p The acceptor doping
     * @param ctype The carrier type (electrons or holes)
     * @return The carrier concentration for thermal bath contacts
     */
    inline double contact_carrier_ohm(double temperature, double doping_n, double doping_p, viennashe::carrier_type_id ctype)
    {
      const double q    = viennashe::physics::constants::q;
      const double kbT  = viennashe::physics::constants::kB * temperature;
      const double phiB = viennashe::physics::built_in_potential(temperature, doping_n, doping_p);

      const double band_edge   = viennashe::physics::get_band_edge(ctype);
      const double band_weight = viennashe::physics::get_band_weight(temperature, ctype);
      const double polarity    = (ctype == viennashe::ELECTRON_TYPE_ID) ? -1.0 : 1.0;

      return band_weight * std::exp( polarity * ( band_edge - q * phiB ) / kbT ) ;
    }


  } //namespace physics
} //namespace viennashe

#endif
