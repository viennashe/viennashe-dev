#ifndef VIENNASHE_MODELS_ALL_BANDGAP_HPP
#define VIENNASHE_MODELS_ALL_BANDGAP_HPP

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

#include <cmath>

#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/math/constants.hpp"
#include "viennashe/models/exception.hpp"

/** @file viennashe/models/bandgap.hpp
    @brief Contains a basic/constant and an extended bandgap model.
 */

namespace viennashe
{
  namespace models
  {

    /**
     * @brief Plain old datatype (POD) of parameters for the bandgap_model_extended
     */
    struct bandgap_model_extended_parameters
    {
      double Eg0;
      double alpha;
      double theta;
      double nu;
      double Tmin;
      double Tmax;

      bandgap_model_extended_parameters() : Eg0(0), alpha(0.0), theta(0.0), nu(0.0), Tmin(0), Tmax(0) { }
    };

    /**
     * @brief Basic bandgap model interface.
     */
    struct bandgap_model
    {
      virtual double operator()(double T) const = 0;
      virtual ~bandgap_model() {}
    };

    /**
     * @brief A simple constant bandgap model. Gets the bandgap and always returns this value.
     */
    class bandgap_model_const : public bandgap_model
    {
    public:

      bandgap_model_const(double Eg) : _Eg(Eg) { }

      double operator()(double) const
      {
        return _Eg;
      }

    private:
      double _Eg;

    };

    /**
     * @brief An elaborate bandgap model, which accounts for the lattice temperature
     */
    class bandgap_model_extended : public bandgap_model
    {
    public:

      bandgap_model_extended(const bandgap_model_extended_parameters & params) : _params(params) { }

      double operator()(double TL) const
      {
        double T = TL;
        if ( T < _params.Tmin ) T = this->_params.Tmin;
        if ( T > _params.Tmax ) T = this->_params.Tmax;

        const double pi = viennashe::math::constants::pi;

        const double prefactor = this->_params.alpha * this->_params.theta * 0.5;
        const double nu = this->_params.nu;
        const double theta = this->_params.theta;

        if ( nu < 0.0 )
          throw viennashe::models::invalid_parameter_exception("bandgap_model_extended: nu must be greater 0.0!", nu);
        if ( nu == 5.0 )
          throw viennashe::models::invalid_parameter_exception("bandgap_model_extended: nu must not be 5.0!", nu);

        const double a1 = (5.0 + nu) / 6.0 * std::pow(pi * 0.5, (2.0 + 0.5 * (nu - 1.0)*(nu - 1.0)));

        const double a2 = (1.0 - nu) * 0.5;

        const double a3 = ((5.0 + nu)*(1.0 + nu)*(1.0 + nu)) / (6.0 * nu + 3.0 * nu * nu);

        const double sum = a1 * std::pow(2.0 * T / theta, 1.0 + nu) + a2 * std::pow(2.0 * T / theta, 2.0 + nu) + a3 * std::pow(2.0 * T / theta, 3.0 + nu);

        const double Eg = _params.Eg0 - prefactor * (std::pow(1.0 + sum + std::pow(2.0 * T / theta, 5.0 + nu), 1.0 / (5.0 + nu)) - 1.0);

        return Eg;
      }

    private:
      bandgap_model_extended_parameters _params;

    };

    /*
    bool test_bandgap_model_extended()
    {
      bandgap_model_extended_parameters params;

      // TEST FOR SILICON
      // Test according to Roland Pässler doi: 10.1002/pssb.200301752
      //   with Data from: W. Bludau, A. Onton, and W. Heinke doi: 10.1063/1.1663501
      params.alpha = viennashe::physics::eV_to_joule(3.22 * 1e-4);
      params.theta = 441.0;
      params.nu = 1.2;
      params.Eg0 = 1.17;

      // TODO: test data needed

      return true;
    }
    */


  } // namespace models

} // namespace viennashe

#endif

