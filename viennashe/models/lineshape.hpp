#ifndef VIENNASHE_MODELS_LINESHAPE_HPP
#define VIENNASHE_MODELS_LINESHAPE_HPP
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


/** @file viennashe/models/lineshape.hpp
    @brief Contains lineshape function "models"
 */

// viennashe
#include "viennashe/forwards.h"

// exceptions
#include "viennashe/models/exception.hpp"

#include "viennashe/math/constants.hpp"
#include "viennashe/physics/constants.hpp"

#include "viennashe/util/checks.hpp"


namespace viennashe
{
  namespace models
  {

    struct lsf_electron_trap_tag {};

    struct lsf_hole_trap_tag {};

    struct lsf_forward_tag {};

    struct lsf_reverse_tag {};

    template < typename TrapTypeTag, typename DirectionTypeTag >
    class lineshape_classic
    {
    public:

      lineshape_classic(double cA, double cB, double qA, double qB, double Et, double epsT2s, double T)
        : cA_(cA), cB_(cB), qA_(qA), qB_(qB), Et_(Et), epsT2s_(epsT2s), T_(T) { }

      void set_temperature(double T) { T_ = T; }

      double operator()(double E) const { return this->operator ()(E, TrapTypeTag(), DirectionTypeTag()); }

    private:

      double cA_, cB_, qA_, qB_, Et_, epsT2s_, T_;

      double operator()(double E, lsf_electron_trap_tag, lsf_forward_tag) const { double ETc = Et_ + epsT2s_; return this->calcLS(cA_, cB_, qB_ - qA_, ETc-E); }
      double operator()(double E, lsf_electron_trap_tag, lsf_reverse_tag) const { double ETc = Et_ + epsT2s_; return this->calcLS(cB_, cA_, qA_ - qB_, E-ETc); }

      double operator()(double E, lsf_hole_trap_tag,     lsf_forward_tag) const { double ETc = Et_ - epsT2s_; return this->calcLS(cA_, cB_, qB_ - qA_, E-ETc); }
      double operator()(double E, lsf_hole_trap_tag,     lsf_reverse_tag) const { double ETc = Et_ - epsT2s_; return this->calcLS(cB_, cA_, qA_ - qB_, ETc-E); }

      double calcLS(double ci, double cf, double qs, double Es) const
      {
        const double pi = viennashe::math::constants::pi;
        const double beta = 1.0/(viennashe::physics::constants::kB * T_);

         double nominator1, nominator2, denominator, q1, q2;

         if (ci == cf)
         {
            q1 = (Es/ci + qs*qs) / (2.*qs);
            return std::sqrt(ci*beta/pi)/2.0 * std::exp(-beta*(ci*q1*q1)) / std::fabs(ci*qs);
         }
         else
         {
            nominator2 = std::sqrt(ci*cf*qs*qs + (ci-cf)*Es);

            if (viennashe::util::is_NaN(nominator2))
              return 0;

            nominator1 = cf*qs;
            denominator = cf-ci;
            q1 = (nominator1 + nominator2)/denominator;
            q2 = (nominator1 - nominator2)/denominator;

            if (std::fabs(q1) < std::fabs(q2))
               return std::sqrt(ci*beta/pi)/2.0 * std::exp(-beta*(ci*q1*q1)) / std::fabs(ci*q1 + cf*(qs-q1));
            else
               return std::sqrt(ci*beta/pi)/2.0 * std::exp(-beta*(ci*q2*q2)) / std::fabs(ci*q2 + cf*(qs-q2));
         }
      }
    };


  } // namespace models
} // namespace viennashe

#endif /* VIENNASHE_MODELS_LINESHAPE_HPP */

