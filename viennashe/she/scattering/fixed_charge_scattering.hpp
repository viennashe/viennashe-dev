#ifndef VIENNASHE_FIXED_CHARGE_SCATTERING_HPP
#define VIENNASHE_FIXED_CHARGE_SCATTERING_HPP

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

#include <limits>

#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/accessors.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/scattering/common.hpp"
#include "viennashe/util/misc.hpp"

/** @file viennashe/she/scattering/fixed_charge_scattering.hpp
    @brief Implements the scattering of carriers with fixed charges.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Trapped charge scattering process
    *
    */
    template <typename DeviceType>
    class fixed_charge_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        explicit fixed_charge_scattering(DeviceType const & device,
                                         viennashe::config const & conf)
        : base_type(device, conf), params_(conf.scattering().fixed_charge())
        { }

        scatter_processes_type operator()(viennagrid_element_id elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(elem, kinetic_energy, ctype);
        }

      protected:

        /** @brief Returns all possible final scattering states for a carrier with initial kinetic energy 'kin_energy' for an edge with vertices v1 and v2
         *
         * @param elem       The element for which to obtain the scattering process
         * @param kinetic_energy  Kinetic energy of the particle before scattering
         * @param ctype      Charge type, either electron or hole
         *
         * @return A vector describing all possible states. Ee use the approximation 'incoming energy equal final energy'.
         */
        scatter_processes_type get(viennagrid_element_id elem,
                                   double kinetic_energy,
                                   viennashe::carrier_type_id ctype) const
        {
          scatter_processes_type result(1);

          // Final energy is always equal to initial energy:
          result[0].initial_energy(kinetic_energy);
          result[0].final_energy(kinetic_energy);
          result[0].rate( getScatteringRate(elem, kinetic_energy, ctype) );

          return result;
        }


        double getScatteringRate(viennagrid_element_id elem,
                                 double kinetic_energy,
                                 viennashe::carrier_type_id ctype) const
        {
           const double temp = base_type::device_.get_lattice_temperature(elem);
           const double NI   = this->get_charged_trap_density(elem, ctype);

           double rate = getScatteringRate(NI, kinetic_energy, temp, ctype);

           // check if the rate is numerically zero
           if ( std::abs(rate) < std::numeric_limits<double>::epsilon() ) rate = 0.0; // assign zero to ensure matrix sparsity

           return rate;
        }


        double get_charged_trap_density(viennagrid_element_id vt,
                                        viennashe::carrier_type_id) const
        {
          const fixed_charge_accessor<DeviceType> fixed_charge(base_type::device_);
          return fixed_charge(vt);
        }


        /*double get_charged_trap_density(viennagrid_element_id facet,
                                        viennashe::carrier_type_id ctype) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<typename base_type::MeshType, FacetType, CellType>::type     CellOnFacetContainer;

          CellOnFacetContainer cells_on_facet(base_type::device_.mesh(), viennagrid::handle(base_type::device_.mesh(), facet));

          CellType const & c1 = *(cells_on_facet.begin());
          CellType const *other_cell_ptr = viennashe::util::get_other_cell_of_facet(base_type::device_.mesh(), facet, c1);

          if (!other_cell_ptr)
            return this->get_charged_trap_density(c1, ctype);

          return std::sqrt(this->get_charged_trap_density(c1, ctype) *
                           this->get_charged_trap_density(*other_cell_ptr, ctype));
        }*/


        double getScatteringRate(double NI, double kinetic_energy, double T, viennashe::carrier_type_id ctype) const
        {
          const double pi   = viennashe::math::constants::pi;
          const double kB   = viennashe::physics::constants::kB;
          const double hbar = viennashe::physics::constants::hbar;
          const double q    = viennashe::physics::constants::q;
          const double eps  = viennashe::materials::si::permittivity();
          double lambda_sq  = 0.0;
          double scattering_rate = 0.0;

          const double norm_k = base_type::conf_.dispersion_relation(ctype).norm_k(kinetic_energy);
          //avoid singularity in the following expression
          if (norm_k <= 0.0 || NI <= 0.0) return 0.0;

          const double prefactor = (2.0 * pi * NI * q * q * q * q) / (hbar * eps * eps);

          lambda_sq = (eps * kB * T) / (q * q * NI) ;

          const double a = 4.0 * lambda_sq * norm_k * norm_k;

          scattering_rate = prefactor * params_.get_fit_factor(ctype)
                            * 0.5 * (std::log(1.0 + a) - (a / (1.0 + a)) ) / 4.0 / std::pow(norm_k, 4.0);

          return scattering_rate;
        }

        scatter_process_id id() const { return FIXED_CHARGE_SCATTERING; }

      private:

        fixed_charge_scattering_parameters params_;
    };

  } //namespace she
} //namespace viennashe


#endif

