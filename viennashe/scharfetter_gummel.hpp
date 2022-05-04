#ifndef VIENNASHE_DD_SCHARFETTER_GUMMEL_HPP
#define VIENNASHE_DD_SCHARFETTER_GUMMEL_HPP

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


#include "viennashe/forwards.h"

#include "viennashe/math/bernoulli.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"

/** @file viennashe/scharfetter_gummel.hpp
    @brief Scharfetter-Gummel stabilization of the continuity equations
*/


namespace viennashe
{
  /** @brief Dispatcher for Scharfetter-Gummel discretization of electron and hold fluxes. */
  class scharfetter_gummel
  {
    public:
      scharfetter_gummel(viennashe::carrier_type_id ctype) : carrier_type_id_(ctype) {}

      /** @brief Scharfetter-Gummel stabilization of the electron and/or hole flux
      *
      * Stabilized way of computing the flux from box i (with electron/hole concentration np_i and potential V_i)
      * to box j (with electron/hole concentration np_j and potential V_j)
      *
      * @param np_i electron/hole concentration in box i
      * @param np_j electron/hole concentration in box j
      * @param V_i potential in box i
      * @param V_j potential in box j
      * @param dx  distance between the two box centers
      * @param mu  The mobility (currently ignored)
      * @param T   temperature
      */
      double operator()(double np_i, double np_j,
                        double V_i, double V_j,
                        double dx, double mu, double T) const
      {
        //double T  = 300.0; // Kelvin
        const double VT = viennashe::physics::get_thermal_potential(T);
        const double D_ij = (V_j-V_i)/VT;

        if (carrier_type_id_ == viennashe::ELECTRON_TYPE_ID)
          return   VT * viennashe::physics::constants::q * mu * VT *        //constant: irrelevant if no recombination is used
                  (np_j*viennashe::math::Bernoulli(D_ij) - np_i*viennashe::math::Bernoulli(-D_ij)) / dx;
        else
          return - VT * viennashe::physics::constants::q * mu * VT *        //constant: irrelevant if no recombination is used
                  (np_j*viennashe::math::Bernoulli(-D_ij) - np_i*viennashe::math::Bernoulli(D_ij)) / dx;
      }

    private:
      viennashe::carrier_type_id carrier_type_id_;
  };



  ////////////////////// Scharfetter-Gummel derivatives with respect to potential




  /** @brief Generic dispatcher for partial derivatives of the Scharfetter-Gummel discretization of electron and hold fluxes. */
  class scharfetter_gummel_dVi
  {
    public:
      scharfetter_gummel_dVi(viennashe::carrier_type_id ctype) : carrier_type_id_(ctype) {}

      /** @brief Scharfetter-Gummel stabilization of the electron flux
      *
      * Stabilized way of computing the flux from box i (with electron/hole concentration np_i and potential V_i)
      * to box j (with electron/hole concentration np_j and potential V_j)
      *
      * @param np_i electron/hole concentration in box i
      * @param np_j electron/hole concentration in box j
      * @param V_i potential in box i
      * @param V_j potential in box j
      * @param dx  distance between the two box centers
      * @param mu  The mobility (currently ignored)
      * @param T   temperature
      */
      double operator()(double np_i, double np_j,
                        double V_i, double V_j,
                        double dx, double mu, double T) const
      {
        const double VT = viennashe::physics::get_thermal_potential(T);
        const double D_ij = (V_j-V_i)/VT;

        if (carrier_type_id_ == viennashe::ELECTRON_TYPE_ID)
          return   viennashe::physics::constants::q * mu * VT * //constant: irrelevant if no recombination is used
                (- np_j*viennashe::math::Bernoulli_dx(D_ij) - np_i*viennashe::math::Bernoulli_dx(-D_ij)) / dx;
        else
          return - viennashe::physics::constants::q * mu * VT * //constant: irrelevant if no recombination is used
                (np_j*viennashe::math::Bernoulli_dx(-D_ij) + np_i*viennashe::math::Bernoulli_dx(D_ij)) / dx;
      }

    private:
      viennashe::carrier_type_id carrier_type_id_;
  };



  /** @brief Generic dispatcher for Scharfetter-Gummel discretization of electron and hold fluxes. */
  class scharfetter_gummel_dVj
  {
    public:
      scharfetter_gummel_dVj(viennashe::carrier_type_id ctype) : carrier_type_id_(ctype) {}

      /** @brief Scharfetter-Gummel stabilization of the electron/hole flux
      *
      * Stabilized way of computing the flux from box i (with electron/hole concentration np_i and potential V_i)
      * to box j (with electron/hole concentration np_j and potential V_j)
      *
      * @param np_i electron concentration in box i
      * @param np_j electron concentration in box j
      * @param V_i potential in box i
      * @param V_j potential in box j
      * @param dx  distance between the two box centers
      * @param mu  The mobility (currently ignored)
      * @param T   temperature
      */
      double operator()(double np_i, double np_j,
                        double V_i, double V_j,
                        double dx, double mu, double T) const
      {
        const double VT = viennashe::physics::get_thermal_potential(T);
        const double D_ij = (V_j-V_i)/VT;

        if (carrier_type_id_ == viennashe::ELECTRON_TYPE_ID)
          return viennashe::physics::constants::q * mu * VT *        //constant: irrelevant if no recombination is used
                (np_j*viennashe::math::Bernoulli_dx(D_ij) + np_i*viennashe::math::Bernoulli_dx(-D_ij)) / dx;
        else
          return viennashe::physics::constants::q * mu * VT *        //constant: irrelevant if no recombination is used
                (np_j*viennashe::math::Bernoulli_dx(-D_ij) + np_i*viennashe::math::Bernoulli_dx(D_ij)) / dx;
      }

    private:
      viennashe::carrier_type_id carrier_type_id_;
  };

} //namespace viennashe

#endif
