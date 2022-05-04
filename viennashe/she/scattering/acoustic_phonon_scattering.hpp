#ifndef VIENNASHE_SHE_SCATTERING_ACOUSTIC_PHONON_HPP
#define VIENNASHE_SHE_SCATTERING_ACOUSTIC_PHONON_HPP
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

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/she/scattering/common.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"


/** @file viennashe/she/scattering/acoustic_phonon_scattering.hpp
    @brief Implements the acoustic phonon scattering processes.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Acoustic phonon scattering process
    *
    * This class provides the parameters for elastic acoustical phonon scattering.
    */
    template <typename DeviceType>
    class acoustic_phonon_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;
        typedef typename base_type::FacetType                   FacetType;
        typedef typename base_type::CellType                    CellType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        explicit acoustic_phonon_scattering(DeviceType const & device,
                                            viennashe::config const & conf) : base_type(device, conf), params_(conf.scattering().acoustic_phonon()) { }

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

        scatter_process_id id() const { return ACOUSTIC_PHONON_SCATTERING; }

      private:

        /** @brief Returns all possible scattering processes involving kinetic energy 'kin_energy' for either an edge or an vertex
         *
         * @param device          A ViennaSHE device object describing the physical device
         * @param elem            Topological element for which to retrieve the scattering process descriptor (either vertex or edge)
         * @param kinetic_energy  Kinetic energy of the particle before scattering
         * @param ctype           Carrier type, either electrons or holes
         *
         * @return A vector describing all possible states. For acoustic phonon scattering we use the approximation 'incoming energy equal final energy'.
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
          result[0].rate( getScatteringRate(base_type::device_.get_lattice_temperature(elem), ctype) );
          result[0].generation_rate(0);

          return result;
        }

        double getScatteringRate(const double T, viennashe::carrier_type_id ctype) const
        {
          const double kB   = viennashe::physics::constants::kB;
          const double hbar = viennashe::physics::constants::hbar;
          const double ul   = params_.get_longitudinal_sound_velocity(ctype);
          const double rho  = params_.get_mass_density(ctype);
          const double E1   = params_.get_deformation_potential(ctype);
          const double pi   = viennashe::math::constants::pi;

          const double scattering_rate = 2.0 * pi * kB * T * E1 * E1 / (hbar * ul * ul * rho);

          if (log_acoustic_phonon_scattering::enabled)
            log::debug<log_acoustic_phonon_scattering>() << "acoustical phonon: " << scattering_rate << std::endl;

          //MUST NOT include density of states and/or N_op
          return params_.get_fit_factor(ctype) * scattering_rate; //_params.get_fit_factor(ctype) * scattering_rate;
        }

        acoustic_phonon_scattering_parameters params_;
    };

  } //namespace she
} //namespace viennashe

#endif

