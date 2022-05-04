#ifndef VIENNASHE_TRAP_LEVEL_HPP
#define VIENNASHE_TRAP_LEVEL_HPP

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

#include "viennashe/physics/physics.hpp"
#include "viennashe/exception.hpp"
#include "log/log.hpp"

/** @file  viennashe/trap_level.hpp
    @brief Contains the definition of a trap level.
*/

namespace viennashe
{

  /** @brief Describes a SRH trap */
  class trap_level
  {
    public:
      explicit trap_level() : collision_cross_section_(0), //inactive by default.
                              density_(0), energy_(0),
                              sign_(-1) { }

      //
      // Collision cross section
      //

      /** @brief Sets the collision cross section (SI units) */
      void collision_cross_section(double ccs)
      {
        if(ccs < 0) throw viennashe::invalid_value_exception("trap_level.collision_cross_section: collision cross sections have to be >= 0 !", ccs);
        collision_cross_section_ = ccs;
      }


      /** @brief Returns the collision cross section (SI units) */
      double collision_cross_section() const { return collision_cross_section_; }

      //
      // Trap density
      //

      /** @brief Returns the trap density for the trap level */
      void density(double d)
      {
        if(d < 0) throw viennashe::invalid_value_exception("trap_level.density: trap densities have to be >= 0 !", d);
        density_ = d;
      }

      /** @brief Returns the trap density for the trap level */
      double density() const { return density_; }
      double charge_sign() const { return sign_; }

      void set_charge_sign(double new_sign) { sign_ = new_sign; }

      void set_donor_like()    { sign_ = -1; }
      void set_acceptor_like() { sign_ = +1; }

      //
      // Trap energy
      //

      /** @brief Returns the trap energy in Joule (zero energy refers to the center of the band gap. */
      void energy(double e) { energy_ = e; }

      /** @brief Returns the trap energy in Joule (zero energy refers to the center of the band gap. */
      double energy() const { return energy_; }

    private:
      double collision_cross_section_;
      double density_;
      double energy_;
      double sign_;
  };


  /** @brief Convenience function for outputting a trap level */
  inline std::ostream & operator<<(std::ostream & os, viennashe::trap_level const & rhs)
  {
    os << "Trap @ " << viennashe::physics::convert::joule_to_eV(rhs.energy()) << " eV "
       << " with css = " << rhs.collision_cross_section() << " m^2 repesentative for "
       << rhs.density() << " m^-3 or m^-2 ";
    return os;
  }

} // namespace viennashe

#endif

