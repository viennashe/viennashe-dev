#ifndef VIENNASHE_PHYSICS_DISPERSION_HPP
#define VIENNASHE_PHYSICS_DISPERSION_HPP

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

#include "viennashe/materials/all.hpp"
#include "viennashe/math/constants.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/materials/dos_RSi.hpp"

#include <math.h>

/** @file viennashe/physics/dispersion.hpp
    @brief Contains the dispersion relations for different materials and different polarities
*/

namespace viennashe
{
  namespace physics
  {

    /** @brief Common interface for band structures
     *
     * Static polymorphism is not used here because this option should be modified from the outside via config files.
     */
    class dispersion_base
    {
      public:
        /** @brief Returns the density of states as a function of kinetic energy (and angles, eventually)
         */
        virtual double density_of_states(double ekin, double theta = 0, double phi = 0) const = 0;

        /** @brief Returns the group velocity as a function of kinetic energy (and angles, eventually)
         */
        virtual double velocity(double ekin, double theta = 0, double phi = 0) const = 0;

        /** @brief Returns the norm of the k-vector as a function of energy  (and angles, eventually). Not possible for all dispersion relations.
         */
        virtual double norm_k(double ekin, double theta = 0, double phi = 0) const = 0;

        /** @brief Returns true if the dispersion relation is isotropic */
        virtual bool is_isotropic() const { return true; }  //be isotropic by default

        /** @brief Clones the dispersion relation. User must ensure that the returned copy is deleted. */
        virtual dispersion_base * clone() const = 0;

        virtual double symmetry_factor() const = 0;

        virtual ~dispersion_base() {}
    };


    /** @brief A proxy object for a dispersion relation. Does NOT take ownership of the provided pointer! */
    class dispersion_proxy
    {
      public:
        dispersion_proxy(dispersion_base const * ptr) : ptr_(ptr) {}

        /** @brief Returns the density of states as a function of kinetic energy (and angles, eventually)
         */
        double density_of_states(double ekin, double theta = 0, double phi = 0) const
        {
          return ptr_->density_of_states(ekin, theta, phi);
        }

        /** @brief Returns the velocity as a function of kinetic energy (and angles, eventually)
         */
        double velocity(double ekin, double theta = 0, double phi = 0) const
        {
          return ptr_->velocity(ekin, theta, phi);
        }

        /** @brief Returns the norm of the k-vector as a function of energy  (and angles, eventually). Not possible for all dispersion relations.
         */
        double norm_k(double ekin, double theta = 0, double phi = 0) const
        {
          return ptr_->norm_k(ekin, theta, phi);
        }

        /** @brief Returns true if the dispersion relation is isotropic */
        bool is_isotropic() const
        {
          return ptr_->is_isotropic();
        }

        dispersion_base const * get() const { return ptr_; }

        double symmetry_factor() const { return ptr_->symmetry_factor(); }

      private:
        dispersion_base const * ptr_;
    };


    /** @brief Parabolic dispersion relation (spherically symmetric)
     *
     * eps = hbar^2 k^2 / (2 m)
     *
     * @tparam MaterialType   A class for the material
     */
    template <typename MaterialType>
    class parabolic_dispersion : public dispersion_base
    {
        typedef parabolic_dispersion   self_type;

      public:
        parabolic_dispersion(viennashe::carrier_type_id ctype) : carrier_type_id_(ctype) {}

        /*double energy(double k_norm) const
        {
          double hbar = constants::hbar;
          double mstar = constants::mstar;
          return hbar*hbar * k_norm * k_norm / (2.0 * mstar);
        }*/

        /** @brief Returns the norm of the k-vector as a function of energy  (and angles, eventually). Not possible for all dispersion relations.
          */
        double norm_k(double energy, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (energy <= 0.0)
            return 0.0;

          double hbar = constants::hbar;
          double mstar = MaterialType::dos_effective_mass(carrier_type_id_);
          return sqrt( 2.0 * energy * mstar / (hbar * hbar) );
        }

        /** @brief Total generalized density of states for one spin direction, considering six-fold symmetry.
        *
        * Since other models may allow for an additional angular dependency, the expression differs by a factor of 4*PI from the conventional expression (4*PI is just the area of the unit sphere)
        *
        * @param eps    The energy (relative to the band edge) for which the generalized DOS should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double density_of_states(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps < 0.0)
            return 0.0;

          const double hbar = constants::hbar;
          const double mstar = MaterialType::dos_effective_mass(carrier_type_id_);
          const double pi = viennashe::math::constants::pi;
          const double prefactor = mstar * sqrt(2.0 * mstar) / (8.0 * pi * pi * pi * hbar * hbar * hbar);

          return prefactor
                * sqrt(eps);
        }

        /** @brief Carrier group velocity
        *
        * @param eps    The energy (relative to the band edge) for which the velocity should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double velocity(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps < 0.0)
            return 0.0;

          double m_c = MaterialType::conductivity_effective_mass(carrier_type_id_);
          double m_d = MaterialType::dos_effective_mass(carrier_type_id_);

          return sqrt(m_c / m_d)  //due to Herring-Vogt transform (see paper by Jin et al.)
                 * sqrt( 2.0 * eps / MaterialType::dos_effective_mass(carrier_type_id_) );
        }

        dispersion_base * clone() const { return new self_type(carrier_type_id_); }

        /** @brief Returns the number of symmetric bands in k-space */
        double symmetry_factor() const
        {
          if (carrier_type_id_ == ELECTRON_TYPE_ID)
            return 6.0;
          return 4.0;
        }

      private:
        viennashe::carrier_type_id carrier_type_id_;
    };


    /** @brief Non-parabolic dispersion relation proposed by the Modena group (spherically symmetric).
     *
     *   eps * (1 + alpha * eps) = hbar^2 * k^2 / (2 * m)
     *
     * @tparam MaterialType   A class for the material
     */
    template <typename MaterialType>
    class modena_dispersion : public dispersion_base
    {
      public:
        modena_dispersion(viennashe::carrier_type_id ctype, double a = 0.5 / viennashe::physics::constants::q) : carrier_type_id_(ctype), alpha(a) {}

        /** @brief Returns the norm of the k-vector as a function of energy  (and angles, eventually). Not possible for all dispersion relations.
         */
        double norm_k(double energy, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (energy < 0.0)
            return 0.0;

          const double hbar = constants::hbar;
          const double mstar = MaterialType::dos_effective_mass(carrier_type_id_);
          return sqrt( 2.0 * energy * (1.0 + alpha * energy) * mstar / (hbar * hbar) );
        }

        /** @brief Generalized density of states
        *
        * Since other models may allow for an additional angular dependency, the expression differs by a factor of 4*PI from the conventional expression (4*PI is just the area of the unit sphere)
        *
        * @param eps    The energy (relative to the band edge) for which the generalized DOS should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double density_of_states(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps > 0.0)
          {
            const double hbar = constants::hbar;
            const double mstar = MaterialType::dos_effective_mass(carrier_type_id_);
            const double pi = viennashe::math::constants::pi;
            const double prefactor = mstar * sqrt(2.0 * mstar) / (8.0 * pi * pi * pi * hbar * hbar * hbar);

            return prefactor * 4.0 * pi
                   * (1.0 + 2.0 * alpha * eps)
                   * sqrt( eps * (1.0 + alpha * eps));
          }

          return 0.0;
        }

        /** @brief Carrier group velocity
        *
        * @param eps    The energy (relative to the band edge) for which the velocity should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double velocity(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps < 0)
            return 0;

          double gamma = eps * (1.0 + alpha * eps);
          double dgam = 1.0 + 2.0 * alpha * eps;  //derivative of gamma w.r.t. eps

          double m_c = MaterialType::conductivity_effective_mass(carrier_type_id_);
          double m_d = MaterialType::dos_effective_mass(carrier_type_id_);

          return sqrt(m_c / m_d)  //due to Herring-Vogt transform (see paper by Jin et al.)
                 * sqrt( 2.0 * gamma / MaterialType::dos_effective_mass(carrier_type_id_)) / dgam;
        }

        dispersion_base * clone() const { return new modena_dispersion(carrier_type_id_, alpha); }

        /** @brief Returns the number of symmetric bands in k-space */
        double symmetry_factor() const
        {
          if (carrier_type_id_ == ELECTRON_TYPE_ID)
            return 6.0;
          return 4.0;
        }

      private:
        viennashe::carrier_type_id carrier_type_id_;
        double alpha;
    };


    /** @brief Non-parabolic isotropic dispersion relation proposed by Jin et al. TED 2011 (uses real DOS)
     *
     *   Uses the real dos for relaxed silicon (dos_RSi.hpp).
     *
     * @tparam MaterialType   A class for the material
     */
    template <typename MaterialType>
    class ext_vecchi_dispersion : public dispersion_base
    {
      public:
        ext_vecchi_dispersion(viennashe::carrier_type_id ctype, double a = 0.5 / viennashe::physics::constants::q) : carrier_type_id_(ctype), alpha(a)  {}

        /** @brief Returns the norm of the k-vector as a function of energy  (and angles, eventually).
         *
         * Uses the dispersion relation from the modena model.
         *
         */
        double norm_k(double energy, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (energy < 0.0) { return 0.0; }

          //throw dispersion_not_invertible_exception("Cannot invert dispersion relation for extended Vecchi model!");
          //return 0;

          double hbar = constants::hbar;
          double mstar = MaterialType::dos_effective_mass(carrier_type_id_);
          return sqrt( 2.0 * energy * mstar / (hbar * hbar) );
        }

        /** @brief Generalized density of states
        *
        * Returns the DOS for relaxed silicon calculated using pseudo potentials
        *
        * @param eps    The energy (relative to the band edge) for which the generalized DOS should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double density_of_states(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps > 0.0)
            return viennashe::materials::dos_relaxed_silicon::get_density_of_states(eps, carrier_type_id_);

          return 0.0;
        }

        /** @brief Carrier group velocity
        *
        * Returns the carrier group velocity for relaxed silicon calculated using pseudo potentials
        *
        * @param eps    The energy (relative to the band edge) for which the velocity should be evaluated
        * @param theta  Inclination
        * @param phi    Azimuth
        */
        double velocity(double eps, double theta = 0, double phi = 0) const
        {
          (void)theta;(void)phi; //prevents unused parameter warnings
          if (eps > 0.0)
          {
            const double m_c = MaterialType::conductivity_effective_mass(carrier_type_id_);
            const double m_d = MaterialType::dos_effective_mass(carrier_type_id_);

            return sqrt(m_c / m_d)  // due to Herring-Vogt transform
                   * viennashe::materials::dos_relaxed_silicon::get_group_velocity(eps, carrier_type_id_);
          }

          return 0.0;
        }

        dispersion_base * clone() const { return new ext_vecchi_dispersion(carrier_type_id_); }

        /** @brief Returns the number of symmetric bands in k-space */
        double symmetry_factor() const
        {
          if (carrier_type_id_ == ELECTRON_TYPE_ID)
            return 6.0;
          return 4.0;
        }

      private:
        viennashe::carrier_type_id carrier_type_id_;
        double alpha;
    };



  } //namespace physics
} //namespace viennashe

#endif
