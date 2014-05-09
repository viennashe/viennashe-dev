#ifndef VIENNASHE_SHE_SCATTERING_OPTICAL_PHONON_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_OPTICAL_PHONON_SCATTERING_HPP

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
#include "viennashe/physics/physics.hpp"
#include "viennashe/she/scattering/common.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"


/** @file viennashe/she/scattering/optical_phonon_scattering.hpp
    @brief Implements the optical phonon scattering processes.
*/

namespace viennashe
{
  namespace she
  {


    /** @brief Optical phonon scattering process
    *
    * This class provides the parameters for inelastic optical phonon scattering.
    * Parameters are a mess at the moment.
    */
    template <typename DeviceType>
    class optical_phonon_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;
        typedef typename base_type::FacetType                   FacetType;
        typedef typename base_type::CellType                    CellType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        /** @brief Constructs the optical phonon scattering object. Takes the phonon energy as parameter
         *
         * @param device               The device
         * @param conf                 The simulator configuration
         * @param energy_grid_spacing  Spacing of the energy grid in order to map phonon energy suitably
         *
         */
        optical_phonon_scattering(DeviceType const & device,
                                  viennashe::config const & conf,
                                  double energy_grid_spacing = 0)
          : base_type(device, conf),
            params_(conf.scattering().optical_phonon()),
            hbaromega_(params_.get_phonon_energy())
        {
          // we have to map the phonon energy to multiples of the grid, otherwise there is an inconsistency with thermal equilibrium
          if (energy_grid_spacing > 0)
          {
            std::size_t grid_steps = static_cast<std::size_t>(hbaromega_ / energy_grid_spacing);
            double error_0 = hbaromega_ - grid_steps * energy_grid_spacing;
            double error_1 = (grid_steps + 1.0) * energy_grid_spacing - hbaromega_;
            if (error_1 < error_0)
              ++grid_steps;

            // Now set new parameters:
            hbaromega_ = grid_steps * energy_grid_spacing;
          }
        }

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

        /////////// scattering to lower energy ///////////


        /** @brief Inelastic energy that is either gained or lost at scattering. */
        double energy() const { return hbaromega_; }

        /** @brief Modifies/Sets the inelastic scattering energy lost or gained
        *
        * Energy is such that initial and final energies lie on grid nodes.
        */
        void energy(double energy)  { hbaromega_ = energy;  }

        scatter_process_id id() const { return OPTICAL_PHONON_SCATTERING; }

      private:

        /** @brief Returns all possible scattering processes involving kinetic energy 'kin_energy' for either an edge or an vertex
         *
         * @param elem            Topological element for which to retrieve the scattering process descriptor (either vertex or edge)
         * @param kinetic_energy  Kinetic energy of the particle before scattering
         * @param ctype           Carrier type, either electrons or holes
         *
         * @return A vector describing all possible states. For optical phonon scattering there are two final scattering states
         */
        template <typename ElementType>
        scatter_processes_type get(ElementType const & elem,
                                   double kinetic_energy,
                                   viennashe::carrier_type_id ctype) const
        {
          scatter_processes_type result(4);

          // Update the optical phonon number
          const double T = base_type::device_.get_lattice_temperature(elem);
          const double N_op = viennashe::physics::get_optical_phonon_number(hbaromega_, T);

          // Emit phonon
          result[0].initial_energy(kinetic_energy);
          result[0].final_energy(kinetic_energy - hbaromega_);
          result[0].rate( getScatteringRate(ctype) * (N_op + 1.0) );
          result[0].generation_rate(0);

          result[1].initial_energy(kinetic_energy + hbaromega_);
          result[1].final_energy(kinetic_energy);
          result[1].rate( getScatteringRate(ctype) * (N_op + 1.0) );
          result[1].generation_rate(0);

          // Capture phonon
          result[2].initial_energy(kinetic_energy);
          result[2].final_energy(kinetic_energy + hbaromega_);
          result[2].rate( getScatteringRate(ctype) * N_op );
          result[2].generation_rate(0);

          result[3].initial_energy(kinetic_energy - hbaromega_);
          result[3].final_energy(kinetic_energy);
          result[3].rate( getScatteringRate(ctype) * N_op );
          result[3].generation_rate(0);

          return result;
        }

        /**
         * @brief Calculates the actual scattering rate
         * @param ctype The carrier type, which does scatter
         * @return The rate in 1/s
         */
        double getScatteringRate(viennashe::carrier_type_id ctype) const
        {
          const double hbar = viennashe::physics::constants::hbar;
          const double pi = viennashe::math::constants::pi;
          const double rho = params_.get_mass_density(ctype);

          const double DtK = params_.get_coupling_constant(ctype);

          const double scattering_rate = pi * DtK * DtK * hbar / (2.0 * rho * hbaromega_);

          //MUST NOT include density of states and/or N_op

          if (log_optical_phonon_scattering::enabled)
            log::debug<log_optical_phonon_scattering>() << "optical phonon: " << scattering_rate << std::endl;

          //return scattering_rate;
          return params_.get_fit_factor(ctype) * scattering_rate;
        }

        optical_phonon_scattering_parameters params_;
        double hbaromega_;
    };


  } //namespace she
} //namespace viennashe


#endif

