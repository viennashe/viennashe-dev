#ifndef VIENNASHE_SHE_SCATTERING_IMPACT_IONIZATION_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_IMPACT_IONIZATION_SCATTERING_HPP

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
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/scattering/common.hpp"

/** @file viennashe/she/scattering/impact_ionization_scattering.hpp
    @brief Implements impact ionization scattering processes.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Impact ionization scattering process
    *
    * This class provides the parameters for impact ionization.
    */
    template <typename DeviceType>
    class impact_ionization_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

        impact_ionization_scattering(DeviceType const & device,
                                     viennashe::config const & conf) : base_type(device, conf), params_(conf.scattering().impact_ionization()) {}

        scatter_processes_type operator()(viennagrid_element_id elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(elem, kinetic_energy, ctype);
        }

        scatter_process_id id() const { return IMPACT_IONIZATION_SCATTERING; }

      private:

        /** @brief Returns all possible scattering processes involving kinetic energy 'kin_energy' for either an edge or an vertex
         *
         * @param device          A ViennaSHE device object describing the physical device
         * @param elem            Topological element for which to retrieve the scattering process descriptor (either vertex or edge)
         * @param kinetic_energy  Kinetic energy of the particle before scattering
         * @param ctype           Carrier type, either electrons or holes
         *
         * @return A vector describing all possible states. For impact ionization we require total energy conservation, in addition there is particle generation.
         */
        template <typename ElementType>
        scatter_processes_type get(ElementType const & /*elem*/,
                                   double kinetic_energy,
                                   viennashe::carrier_type_id ctype) const
        {
          scatter_processes_type result(2);

          // What is II?
          //  An electron (primary) in the CB scatters with an electron (secondary) in the VB.
          //  The primary electron looses energy (from eps to (eps-e_gap)/3 and stays in the CB.
          //  The secondary electron is lifted from the VB over the bandgap into the CB: energy = (eps-e_gap)/3 (same as primary!).
          //  When the secondary electron is lifted to the CB a hole is left behind in the VB.
          //  This is the secondary hole. Its energy is (eps-e_gap)/3 (same as primary!)
          //  Of course we'd need also to account for the momentum (k-vector) conservation, but it turned out that the
          //   (eps-e_gap)/3 rule is a good approximation [Jungemann et al.].
          //  Thus we need to generate both carrier types and scatter the primary carrier!

          // out-scattering:
          result[0].initial_energy( kinetic_energy );
          result[0].final_energy( (kinetic_energy - viennashe::materials::si::band_gap()) / 3.0 );
          result[0].rate( getScatteringRate(kinetic_energy,
                                            ctype) );
          result[0].generation_rate(0.0);

          // in-scattering (generates new carriers):
          result[1].initial_energy( 3.0 * kinetic_energy + viennashe::materials::si::band_gap() );
          result[1].final_energy( kinetic_energy );
          result[1].rate( getScatteringRate(3.0 * kinetic_energy + viennashe::materials::si::band_gap(),
                                            ctype) );
          result[1].generation_rate( get_secondary_carrier_generation_rate(kinetic_energy) );

          return result;
        }


        double get_secondary_carrier_generation_rate(const double kinetic_energy) const
        {
          // This makes the approximation apparent ... after scattering all three carriers have the same energy.

          return   approximated_scattering_factor_n(kinetic_energy) // from electron => electron + hole
                 + approximated_scattering_factor_p(kinetic_energy); // from hole => hole + electron
        }

        /** @brief Approximation of II according to Jungemann and Meinerzhagen */
        double approximated_scattering_factor_n(const double kinetic_energy) const
        {
          const double a_low  = 1.49e11; // 1/s
          const double a_high = 1.13e12; // 1/s
          const double b_low  = 1.128; // eV
          const double b_high = 1.572; // eV

          const double kinetic_energy_eV = kinetic_energy / viennashe::physics::constants::q; // Joule => eV

          double scat = 0.0;

          if ( 1.128 <= kinetic_energy_eV )
          {
            log::debug<log_impact_ionization_scattering>() << "kinetic_energy = " << kinetic_energy_eV << " eV " << std::endl;

            if ( kinetic_energy_eV < 1.750)
            {
              scat = a_low  * ( kinetic_energy_eV - b_low ) * ( kinetic_energy_eV - b_low ) * ( kinetic_energy_eV - b_low ) ;
            }
            else
            {
              scat = a_high * ( kinetic_energy_eV - b_high ) * ( kinetic_energy_eV - b_high ) ;
            }
          }
          return scat;
        }

        /** @brief Impact ionization for holes  */
        double approximated_scattering_factor_p(const double kinetic_energy) const
        {
          double scat = 0;

          const double kinetic_energy_eV = kinetic_energy / viennashe::physics::constants::q; // Joule => eV

          // Fischetti
          const double eth1 = 1.125; // eV
          const double eth2 = 1.45;  // eV
          const double A = 0.002e12; // 1/s
          const double B = 1.0e12; // 1/s

          // from MONJU
          //const double eth1 = 1.125; // eV
          //const double eth2 = 1.45;  // eV
          //const double A = 0.5970318495776689; // 1/s // CHECK THIS
          //const double B = 446661.31082650006; // 1/s // CHECK THIS

          // Note that in the following the energy is in eV and everything in std::pow is normed to 1 eV
          if ( kinetic_energy_eV > eth1)
          {
            scat  = A * std::pow( (kinetic_energy_eV - eth1) /* / 1 eV */, 6 );
          }
          if ( kinetic_energy_eV > eth2)
          {
            scat += B * std::pow( (kinetic_energy_eV - eth2) /* / 1 eV */, 4);
          }

          return scat;
        }


        double getScatteringRate(const double kinetic_energy,
                                 viennashe::carrier_type_id ctype) const
        {
          const double eps = kinetic_energy;

          const double summed_dos = base_type::conf_.dispersion_relation(ctype).density_of_states(eps);// TODO: Sum DOS over ALL bands !!

          //avoid singularity in the following expression
          if (summed_dos <= 0.0)
          {
            //log::warn<log_impact_ionization_scattering>() << " Summed DOS <= 0.0 ## eps = " << eps << std::endl;
            return 0.0;
          }

          if (ctype == ELECTRON_TYPE_ID)
          {
            log::debug<log_impact_ionization_scattering>() << " approximated_scattering_factor(kinetic_energy) / summed_dos = " << approximated_scattering_factor_n(kinetic_energy) / summed_dos << std::endl;
            return approximated_scattering_factor_n(kinetic_energy) / summed_dos;
          }

          log::debug<log_impact_ionization_scattering>() << " approximated_scattering_factor(kinetic_energy) / summed_dos = " << approximated_scattering_factor_p(kinetic_energy) / summed_dos << std::endl;
          return approximated_scattering_factor_p(kinetic_energy) / summed_dos;
        }

        impact_ionization_scattering_parameters params_;
    };

  } //namespace she
} //namespace viennashe

#endif

