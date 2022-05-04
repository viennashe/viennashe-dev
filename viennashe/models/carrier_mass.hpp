#ifndef VIENNASHE_MODELS_ALL_CARRIER_MASS_HPP
#define VIENNASHE_MODELS_ALL_CARRIER_MASS_HPP

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
#include "viennashe/exception.hpp"

/** @file viennashe/models/carrier_mass.hpp
    @brief Contains the basic carrier mass models (Taylor series).
 */

namespace viennashe
{
  namespace models
  {

    /** @brief POD with parameters for the simple carrier mass model (constant value model) */
    struct carrier_mass_simple_model_parameters
    {
      // The mass of the electron
      double mass_electron;
      // The mass of the hole
      double mass_hole;

      carrier_mass_simple_model_parameters() : mass_electron(0), mass_hole(0) { }
    };

    /** @brief POD with parameters for the full model */
    struct carrier_mass_full_model_parameters
    {
      // Taylor series (with lattice temperature) coefficients for the electron mass
      double mass_electron_offset, mass_electron_linear, mass_electron_quadradic;
      // Taylor series (with lattice temperature) coefficients for the hole mass
      double mass_hole_offset, mass_hole_linear, mass_hole_quadradic;

      carrier_mass_full_model_parameters() : mass_electron_offset(0), mass_electron_linear(0), mass_electron_quadradic(0),
      mass_hole_offset(0), mass_hole_linear(0), mass_hole_quadradic(0) { }
    };


    /** @brief The carrier mass model interface */
    struct carrier_mass_model
    {
      /**
       * @brief The functor interface
       * @param TL The lattice temperature
       * @param ctype The carrier type
       */
      virtual double operator()(double TL, viennashe::carrier_type_id ctype) const = 0;
      virtual ~carrier_mass_model() {}
    }; // carrier_mass_model


    /** @brief The simple (constant) carrier mass model */
    class carrier_mass_simple_model : public carrier_mass_model
    {
    protected:
      const carrier_mass_simple_model_parameters _params;

    public:
      carrier_mass_simple_model(const carrier_mass_simple_model_parameters & params) : _params(params) { }

      double operator()(double TL, viennashe::carrier_type_id ctype) const
      {
        (void)TL; //prevent unused parameter warnings
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return _params.mass_electron;
        else if (ctype == viennashe::HOLE_TYPE_ID)
          return _params.mass_hole;
        else
          throw viennashe::carrier_type_not_supported_exception("carrier_mass_simple_model only supports holes and electrons");
      }

    }; // carrier_mass_simple_model


    /** @brief The full carrier mass model, which uses a Taylor series over lattice temperature */
    class carrier_mass_full_model : public carrier_mass_model
    {
    protected:
      const carrier_mass_full_model_parameters _params;

    public:
      carrier_mass_full_model(const carrier_mass_full_model_parameters & params)
      : _params(params) { }

      double operator()(double TL, viennashe::carrier_type_id ctype) const
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return this->for_electrons(TL);
        else if (ctype == viennashe::HOLE_TYPE_ID)
          return this->for_holes(TL);
        else
          throw viennashe::carrier_type_not_supported_exception("carrier_mass_simple_model only supports holes and electrons");
      }

    private:

      double for_electrons(double TL) const
      {
        const double TL300 = TL / 300.0;
        return _params.mass_electron_offset
            + _params.mass_electron_linear * TL300;
      }

      double for_holes(double TL) const
      {
        const double TL300 = TL / 300.0;
        return _params.mass_hole_offset
            + _params.mass_hole_linear * TL300
            + _params.mass_hole_quadradic * TL300 * TL300;
      }

    }; // carrier_mass_full_model



  } // namespace models

} // namespace viennashe

#endif

