#ifndef VIENNASHE_SHE_SCATTERING_IMPURITY_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_IMPURITY_SCATTERING_HPP

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
#include "viennashe/math/constants.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/scattering/common.hpp"

/** @file viennashe/she/scattering/impurity_scattering.hpp
    @brief Implements the ionized impurity scattering processes
*/

namespace viennashe
{
  namespace she
  {


    /** @brief Ionized impurity scattering process
    *
    * This class provides the parameters for ionized impurity scattering.
    */
    template <typename DeviceType>
    class ionized_impurity_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        explicit ionized_impurity_scattering(DeviceType const & device,
                                             viennashe::config const & conf) : base_type(device, conf), params_(conf.scattering().ionized_impurity()) {  }

        scatter_processes_type operator()(viennagrid_element_id elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return this->get(elem, kinetic_energy, ctype);
        }

        scatter_process_id id() const { return IMPURITY_SCATTERING; }

      private:

        /** @brief Returns all possible scattering processes involving kinetic energy 'kin_energy' for either an edge or an vertex
         *
         * @param elem            Topological element for which to retrieve the scattering process descriptor (either vertex or edge)
         * @param kinetic_energy  Kinetic energy of the particle before scattering
         * @param ctype           Carrier type, either electrons or holes
         *
         * @return A vector describing all possible states. For ionized_impurity_scattering we use the approximation 'incoming energy equal final energy'.
         */
        template <typename ElementType>
        scatter_processes_type get(ElementType const & elem,
                                   double kinetic_energy,
                                   viennashe::carrier_type_id ctype) const
        {
          scatter_processes_type result(1);

          // Final energy is always equal to initial energy:
          result[0].initial_energy(kinetic_energy);
          result[0].final_energy(kinetic_energy);
          result[0].rate( getScatteringRate(elem, kinetic_energy, ctype) );
          result[0].generation_rate(0);

          return result;
        }


        template <typename ElementType>
        double getScatteringRate(ElementType const & elem,
                                 double kinetic_energy,
                                 viennashe::carrier_type_id ctype) const
        {
           double ND = base_type::device_.get_doping_n(elem);
           double NA = base_type::device_.get_doping_p(elem);
           const double temp = base_type::device_.get_lattice_temperature(elem);

           double NI = NA + ND;

           return getScatteringRate(NI, kinetic_energy, temp, ctype);
        }

        double getScatteringRate(double NI,
                                 double kinetic_energy, double T, viennashe::carrier_type_id ctype) const
        {
          const double kB = viennashe::physics::constants::kB;
          const double hbar = viennashe::physics::constants::hbar;
          const double q = viennashe::physics::constants::q;
          const double eps = viennashe::materials::si::permittivity();
          const double pi = viennashe::math::constants::pi;

          const double prefactor = (2.0 * pi * NI * q * q * q * q) / (hbar * eps * eps);

          const double norm_k = base_type::conf_.dispersion_relation(ctype).norm_k(kinetic_energy);

          //avoid singularity in the following expression
          if (norm_k <= 0.0)
            return 0.0;

          const double lambda_sq = (eps * kB * T) / (q * q * NI );
          const double a = 4.0 * lambda_sq * norm_k * norm_k;

          const double scattering_rate = prefactor
                                         * 0.5 * (std::log(1.0 + a) - (a / (1.0 + a)) ) / 4.0 / std::pow(norm_k, 4.0);

          if (scattering_rate < 0.0) //the above formula may in some corner cases become negative, so we truncate it
            return 0.0;

          return scattering_rate;
        }

        ionized_impurity_scattering_parameters  params_;

    };

  } //namespace she
} //namespace viennashe

#endif

