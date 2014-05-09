#ifndef VIENNASHE_MODELS_TUNNELING_HPP
#define VIENNASHE_MODELS_TUNNELING_HPP
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

/** @file viennashe/models/tunneling.hpp
    @brief Contains tunneling models
 */

#include <cmath>

// viennashe
#include "viennashe/forwards.h"

// exceptions
#include "viennashe/models/exception.hpp"

#include "viennashe/physics/constants.hpp"


namespace viennashe
{
  namespace models
  {

    namespace detail
    {
      double wkb_tunneling_impl(double E, double psi_i, double psi_t, double x, double F,
                                double mtun = 0.33)
      {
        const double m0   = viennashe::physics::constants::mass_electron;
        const double hbar = viennashe::physics::constants::hbar;
        const double q0   = viennashe::physics::constants::q;

        const double C  = 4 * std::sqrt(2*m0)/(3*hbar*q0);
        const double Cs = 2 * std::sqrt(2*m0)/hbar;

        if (std::abs(F) < 5e6)
        {
          if (psi_i > E)
            return std::exp(-Cs*std::sqrt(mtun*(psi_i-E))*x);
          return 1.0;
        }

        if (psi_i >= E && psi_t >= E)
        {
          return std::exp( -C * std::sqrt(mtun)/F * (std::pow((psi_t - E), 1.5) - std::pow((psi_i - E), 1.5)) );
        }
        else if (psi_t >= E && psi_i <= E)
        {
          return std::exp( -C * std::sqrt(mtun)/F * (std::pow((psi_t - E), 1.5)));
        }
        else if (psi_i >= E && psi_t <= E)
        {
          return std::exp(  C * std::sqrt(mtun)/F * (std::pow((psi_i - E), 1.5)));
        }

        return 1.0;
      }
    } // namespace detail

    /** @brief A functor interface to the WKB tunneling coefficient for oxide barriers */
    class wkb_oxide_barrier_tunneling
    {
    public:

      wkb_oxide_barrier_tunneling(double psi_i, double psi_t, double x, double F = 0.0)
        : psi_i_(psi_i), psi_t_(psi_t), x_(x), F_(F) {}

      void set_x(double x) { x_ = x; }
      void set_F(double F) { F_ = F; }

      /**
       * @brief The main user function to get the WKB coefficient |T|^2 ( 0 <= |T|^2 <= 1).
       * @param E  The energy of the charge carrier to tunnel (Joule)
       * @return The WKB coefficient |T|^2 ( 0 <= |T|^2 <= 1)
       */
      double operator()(double E) const
      {
        return viennashe::models::detail::wkb_tunneling_impl(E, psi_i_, psi_t_, x_, F_);
      }

    private:
      double psi_i_; // Hight of the barrier in Joule from the band edge on the far side of the barrier wrspt. the trap
      double psi_t_; // Hight of the barrier in Joule from the band edge on the side of the trap
      double x_; // The thickness of the tunneling barrier (meter)
      double F_; // The electric field at the (most probable) origin of the carrier (V/m)
    };


  } // namespace models
} // namespace viennashe

#endif /* VIENNASHE_MODELS_TUNNELING_HPP */

