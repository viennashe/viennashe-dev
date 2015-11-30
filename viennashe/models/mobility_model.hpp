#ifndef VIENNASHE_MODELS_MOBILITY_MODEL_HPP
#define VIENNASHE_MODELS_MOBILITY_MODEL_HPP

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
#include <stdexcept>
#include <utility>

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/models/mobility_parameters.hpp"
#include "viennashe/models/exception.hpp"

// viennagrid
#include "viennagrid/viennagrid.h"


/** @file viennashe/models/mobility_model.hpp
    @brief Contains basic mobility models for lattice, impurity, field and surface mobility changes
 */


namespace viennashe
{
  namespace models
  {

    /** @brief This namespace contains models, which are only for the drift diffusion transport model */
    namespace dd
    {

      namespace mobility_detail
      {

        class mobility_lattice_scattering
        {
        public:

          mobility_lattice_scattering(const lattice_scattering & params) : _params(params) { }

          double operator()(double mu, double T) const
          {
            if ( _params.enabled )
            {
              return this->calculateMobilityLatticeScattering(mu, T);
            }
            else
            {
              return mu;
            }
          }

        private:
          lattice_scattering _params;

          double calculateMobilityLatticeScattering(double mu, double T) const
          {
            if ( _params.T_ref == 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.T_ref = 0 and therefore invalid");

            return mu * pow(T / _params.T_ref, -_params.alpha);
          }
        }; // mobility_lattice_scattering

        class mobility_impurity_scattering
        {
        public:

          mobility_impurity_scattering(const impurity_scattering & params) : _params(params) { }

          double operator()(double mu, double total_doping) const
          {
            if ( _params.enabled )
            {
              return this->calculateMobilityIonizedImpurityScattering(mu, total_doping);
            }
            else
            {
              return mu;
            }
          }

        private:
          impurity_scattering _params;

          double calculateMobilityIonizedImpurityScattering(double mu, double total_doping) const
          {
            if ( _params.N_ref == 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.N_ref = 0 and therefore invalid");
            if ( _params.mu_min < 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.mu_min < 0 and therefore invalid");
            if ( total_doping <= 0.0 )
              throw viennashe::models::invalid_parameter_exception("total_doping <= 0 and therefore invalid");

            double N_sum = 1.0 * total_doping;

            // Caughey and Thomas
            double N_adjust = pow(N_sum / _params.N_ref, _params.alpha);
            return _params.mu_min + (mu - _params.mu_min) / (1 + N_adjust);
          }
        }; // mobility_impurity_scattering

        class mobility_surface_scattering
        {
        public:

          mobility_surface_scattering(const surface_scattering & params) : _params(params) { }

          double operator()(double mu, double distance_to_insulator, double E_n) const
          {
            if ( _params.enabled )
            {
              return this->calculateMobilitySurfaceScattering(mu, distance_to_insulator, E_n);
            }
            else
            {
              return mu;
            }
          }

        private:
          surface_scattering _params;

          double calculateMobilitySurfaceScattering(double mu, double distance_to_insulator, double E_n) const
          {
            if ( _params.depth_ref == 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.depth_ref = 0 and therefore invalid");
            if ( _params.E_ref == 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.E_ref = 0 and therefore invalid");
            if ( _params.mu_ref < 0.0 )
              throw viennashe::models::invalid_parameter_exception("params.mu_ref < 0 and therefore invalid");
            if ( distance_to_insulator < 0.0 )
              throw viennashe::models::invalid_parameter_exception("distance_to_next_insulator < 0 and therefore invalid");

            double depth = distance_to_insulator;

            // Selberherr-Formula
            double depth_dependence = 2.0 * exp(-pow(depth / _params.depth_ref, 2)) / (1 + exp(-2.0 * pow(depth / _params.depth_ref, 2)));
            double field_dependence = pow(E_n / _params.E_ref, _params.gamma_ref);
            return ( _params.mu_ref + (mu - _params.mu_ref) * (1.0 - depth_dependence)) / (1.0 + depth_dependence * field_dependence);
          }

        }; // mobility_surface_scattering

        class mobility_field_dependence
        {
        public:

          mobility_field_dependence(const field_dependence & params) : _params(params) { }

          double operator()(double mu, double T, double E_n) const
          {
            if ( _params.enabled )
            {
              return this->calculateFieldDependentMobility(mu, T, E_n);
            }
            else
            {
              return mu;
            }
          }

        private:
          field_dependence _params;

          double get_saturation_velocity(double T) const
          {
            return _params.vsat300 / (1.0 + _params.vsat300C * (T / 300.0 - 1.0));
          }

          double calculateFieldDependentMobility(double mu, double T, double E_n) const
          {
            double F = E_n;
            double v_sat = this->get_saturation_velocity(T);
            double beta = _params.beta;

            if ( v_sat <= 0.0 )
              throw viennashe::models::invalid_parameter_exception("v_sat <= 0 and therefore invalid");
            if ( T <= 0.0 )
              throw viennashe::models::invalid_parameter_exception("T <= 0 and therefore invalid");

            double h = pow(pow(2.0 * F * mu / v_sat, beta) + 1.0, beta) + 1.0;

            return 2.0 * mu / h;
          }


        }; // mobility_field_dependence

      } // namespace mobility_detail


      /** @brief The main mobility model. Contains submodels for lattice, impurity, field and
       *         surface scattering related mobility changes. For drift diffusion only! */
      template < typename DeviceType >
      class mobility
      {
      public:
        typedef double value_type;

        mobility(const DeviceType & device, const mobility_paramters & params) : _params(params), _device(device) { }

        mobility(const mobility & mob) : _params(mob._params), _device(mob._device) { }

        /**
         * @brief Returns the mobility along a facet
         * @param facet The facet
         * @param potential An accessor (for cells) to the electrostatic potential
         * @return A mobility in SI units
         */
        /*template < typename PotentialAccessor >
        value_type operator()(const FacetType & facet, PotentialAccessor const & potential) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
          typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

          CellOnFacetContainer cells_on_facet(_device.mesh(), viennagrid::handle(_device.mesh(), facet));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          CellType const & c1 = *cofit;
          ++cofit;
          CellType const & c2 = *cofit;

          return this->operator ()(c1, c2, potential);
        }*/

        /**
         * @brief Returns the mobility from a vertex to another vertex. Usage of the other, edge related, operator() is recommended.
         * @param c1 The source vertex
         * @param c2 The sink vertex
         * @param potential An accessor (for vertices) to the electrostatic potential
         * @return A mobility in SI units
         */
        template < typename PotentialAccessor >
        value_type operator()(viennagrid_element_id c1, viennagrid_element_id c2, PotentialAccessor const & potential) const
        {
          mobility_detail::mobility_lattice_scattering  lattice(_params.lattice);
          mobility_detail::mobility_impurity_scattering impurity(_params.impurity);
          //mobility_detail::mobility_surface_scattering  surface(_params.surface);
          mobility_detail::mobility_field_dependence    field(_params.field);

          const double TL = _device.get_lattice_temperature(c1);
          const double total_doping_on_cell = (_device.get_doping_n(c1) + _device.get_doping_p(c1));

          std::vector<double> centroid_1(3), centroid_2(3);
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(_device.mesh(), c1, &(centroid_1[0])));
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(_device.mesh(), c2, &(centroid_2[0])));
          double edge_len;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(3, &(centroid_1[0]), &(centroid_2[0]), &edge_len));

          const double    Emag   = ( -(potential(c2) - potential(c1)) / edge_len);

          double mu = _params.mu0;
          if (mu <= 0.0) throw viennashe::models::invalid_parameter_exception("The constant mobility parameter 'mu0' is smaller or equal 0!");

          mu = lattice (mu, TL);
          mu = impurity(mu, total_doping_on_cell );

          /*if ( (_params.surface.enabled ) && _device.has_distance_to_next_insulator(vt))
          {
            const PointType n      = ( viennagrid::centroid(c2) - viennagrid::centroid(c1) ) / edge_len;
            const PointType E_edge = Emag * n;
            const double distance         = _device.distance_to_next_insulator(vt);
            const PointType dist_vector_n = _device.vector_to_next_insulator(vt) / distance; // normalized distance vector

            const double En = viennagrid::inner_prod(E_edge, dist_vector_n); // projection of E on distance to insulator vector

            mu = surface (mu, distance, En);
          }*/

          mu = field (mu, TL, Emag);

          return mu;
        }

      private:
        const mobility_paramters  _params;
        const DeviceType & _device;

      }; // mobility

    } //namespace dd

  } // namespace models

} //namespace viennashe

#endif
