#ifndef VIENNASHE_MODELS_MOBILITY_HPP
#define VIENNASHE_MODELS_MOBILITY_HPP

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

// viennashe
#include "viennashe/forwards.h"

#include "viennashe/models/mobility_parameters.hpp"
#include "viennashe/models/mobility_model.hpp"

namespace viennashe
{
  namespace models
  {
    /** @brief Compiletime evaluation namespace */
    namespace result_of
    {
      /** @brief Compiletime mobility-type getter */
      template < typename DeviceType >
      struct mobility_type
      {
        typedef typename viennashe::models::dd::mobility<DeviceType> type;
      };

    } // result_of

    /**
     * @brief Creates a new mobility model using the given parameters
     * @param device The device
     * @param params The mobility model parameters
     * @return A viennashe::models::dd::mobility model
     */
    template < typename DeviceType >
    viennashe::models::dd::mobility<DeviceType>
      create_mobility_model(DeviceType const & device,
                            viennashe::models::dd::mobility_paramters const & params)
    {
      return viennashe::models::dd::mobility<DeviceType>(device, params);
    }

    /**
     * @brief Returns a mobility model, which always yields the same mobility
     * @param device
     * @param mu The mobility you like the model to yield
     * @return A viennashe::models::dd::mobility model, where all submodels have been disabled
     */
    template < typename DeviceType >
    viennashe::models::dd::mobility<DeviceType> create_constant_mobility_model(DeviceType const & device, double mu)
    {
      viennashe::models::dd::mobility_paramters params;
      params.mu0 = mu;
      params.field.enabled    = false;
      params.lattice.enabled  = false;
      params.impurity.enabled = false;
      params.surface.enabled  = false;

      return viennashe::models::dd::mobility<DeviceType>(device, params);
    }


  } // models
} // viennashe

#endif /* VIENNASHE_MODELS_MOBILITY_HPP */

