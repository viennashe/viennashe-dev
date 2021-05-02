#ifndef VIENNASHE_MODELS_DD_MOBILITY_PARAMETERS_H
#define VIENNASHE_MODELS_DD_MOBILITY_PARAMETERS_H
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

/** @file viennashe/models/mobility_parameters.hpp
    @brief Contains PODs for the mobility parameters of the mobility models
 */


namespace viennashe
{
  namespace models
  {
    namespace dd
    {
      /** @brief Hides implementation details */
      namespace mobility_detail
      {
        struct lattice_scattering
        {
          double alpha;
          double T_ref;
          bool enabled;

          lattice_scattering() : alpha(0), T_ref(0), enabled(false) { }
        };

        struct impurity_scattering
        {
          double mu_min;
          double alpha;
          double N_ref;
          bool enabled;

          impurity_scattering() : mu_min(0), alpha(0), N_ref(0), enabled(false) { }
        };

        struct surface_scattering
        {
          double mu_ref;
          double E_ref;
          double depth_ref;
          double gamma_ref;
          bool enabled;

          surface_scattering() : mu_ref(0), E_ref(0), depth_ref(0), gamma_ref(0), enabled(false) { }
        };

        struct field_dependence
        {
          double beta;
          double vsat300;
          double vsat300C;
          bool enabled;

          field_dependence() : beta(0), vsat300(0), vsat300C(0), enabled(false) { }
        };
      } // mobility_detail

      /** @brief The combined POD for the mobility parameters */
      struct mobility_paramters
      {
        double mu0;
        mobility_detail::lattice_scattering  lattice;
        mobility_detail::impurity_scattering impurity;
        mobility_detail::surface_scattering  surface;
        mobility_detail::field_dependence    field;

        mobility_paramters() : mu0(0.1430) { }
      };


    } // namespace models
  } // namespace dd

} // namespace viennashe

#endif

