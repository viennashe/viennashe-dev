#ifndef VIENNASHE_TRAPPED_CHARGE_SCATTERING_HPP
#define VIENNASHE_TRAPPED_CHARGE_SCATTERING_HPP

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

#include <limits>

// ViennaGrid:

#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/mesh/coboundary_iteration.hpp"

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/scattering/common.hpp"
#include "viennashe/util/misc.hpp"


/** @file viennashe/she/scattering/trapped_charge_scattering.hpp
    @brief Implements the trapped charge scattering processes.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Trapped charge scattering process
    *
    *
    */
    template <typename DeviceType, typename TimeStepQuantitiesT>
    class trapped_charge_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;
        typedef typename base_type::FacetType                   FacetType;
        typedef typename base_type::CellType                    CellType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        explicit trapped_charge_scattering(DeviceType const & device,
                                           viennashe::config const & conf,
                                           TimeStepQuantitiesT const & quantities)
          : base_type(device, conf), params_(conf.scattering().trapped_charge()), quantities_(quantities)
        { }

        scatter_processes_type operator()(CellType const & elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(elem, kinetic_energy, ctype);
        }

        scatter_processes_type operator()(FacetType const & elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(elem, kinetic_energy, ctype);
        }

        scatter_process_id id() const { return TRAPPED_CHARGE_SCATTERING; }

      protected:

        /** @brief Returns all possible final scattering states for a carrier with initial kinetic energy 'kin_energy' for an edge with vertices v1 and v2
         *
         * @param elem The element for which to obtain the scattering process
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

          return result;
        }


        template <typename ElementType>
        double getScatteringRate(ElementType const & elem,
                                 double kinetic_energy,
                                 viennashe::carrier_type_id ctype) const
        {
           const double temp = base_type::device_.get_lattice_temperature(elem);

           const double NI = std::abs(this->get_charged_trap_density(elem, ctype));

           double rate = get_scattering_rate(NI, kinetic_energy, temp, ctype);

           // check if the rate is numerically zero
           if ( std::abs(rate) < std::numeric_limits<double>::epsilon() ) rate = 0.0; // assign zero to ensure matrix sparsity

           return rate;
        }


        double get_charged_trap_density(CellType const & cell,
                                        viennashe::carrier_type_id) const
        {
          typedef typename DeviceType::trap_level_container_type     trap_level_container_type;
          typedef typename trap_level_container_type::const_iterator trap_iterator_type;

          double N = 0.0;

          trap_level_container_type const & traps = this->device_.get_trap_levels(cell);

          const std::size_t num_trap_unknowns = quantities_.num_trap_unknown_indices(cell);

          if (num_trap_unknowns <= 0)
            return 0.0;

          if (num_trap_unknowns != traps.size())
            throw viennashe::invalid_value_exception("The number of traps configured in the device does not match the number of unknowns for traps!", static_cast<double>(num_trap_unknowns));

          std::size_t index = 0;
          for (trap_iterator_type tit = traps.begin(); tit != traps.end(); ++tit, ++index)
          {
            const double occupancy = quantities_.trap_occupancy(cell, index);
            N +=  tit->charge_sign() * occupancy * tit->density() ;
          } // for trap levels

          return N;

        }

        double get_charged_trap_density(FacetType const & facet,
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
        }


        double get_scattering_rate(double NI, double kinetic_energy, double T, viennashe::carrier_type_id ctype) const
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

          //if(scattering_rate != 0) log::debug() << "trapped charge scat    " << kinetic_energy << "  " << scattering_rate << std::endl;

          return scattering_rate;
        }

      private:

        trapped_charge_scattering_parameters params_;

        TimeStepQuantitiesT const & quantities_;

    };

  } //namespace she
} //namespace viennashe


#endif /* VIENNASHE_TRAPPED_CHARGE_SCATTERING_HPP */

