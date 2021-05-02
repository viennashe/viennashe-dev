#ifndef VIENNASHE_SHE_SCATTERING_SURFACE_ACOUSTIC_PHONON_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_SURFACE_ACOUSTIC_PHONON_SCATTERING_HPP

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

#ifdef VIENNASHE_USE_DISABLED_CODE

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

// viennagrid
#include "viennagrid/forwards.h"
#include "viennagrid/algorithm/inner_prod.hpp"
#include "viennagrid/algorithm/norm.hpp"

/** @file viennashe/she/surface_acoustic_phonon_scattering.hpp
    @brief Implements the surface acoustic phonon scattering processes.
 */

namespace viennashe
{
  namespace she
  {

    /** @brief Trapped charge scattering process
     *
     *
     */
    template <typename DeviceType, typename ControllerType, typename ElectricFieldAccessor>
    class surface_acoustic_phonon_scattering : public scattering_base<DeviceType, ControllerType>
    {
        typedef scattering_base<DeviceType, ControllerType>     base_type;
        typedef typename base_type::FacetType                   FacetType;
        typedef typename base_type::CellType                    CellType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

      explicit surface_acoustic_phonon_scattering(const viennashe::materials::surface_acoustic_phonon_scattering_parameters & params,
          ElectricFieldAccessor const & Efield)
       : _params(params), _Efield(Efield) { }


        scatter_processes_type operator()(DeviceType const & device,
                                          ControllerType const & controller,
                                          VertexType const & elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(device, controller, elem, kinetic_energy, ctype);
        }

        scatter_processes_type operator()(DeviceType const & device,
                                          ControllerType const & controller,
                                          EdgeType const & elem,
                                          double kinetic_energy,
                                          viennashe::carrier_type_id ctype) const
        {
          return get(device, controller, elem, kinetic_energy, ctype);
        }

        scatter_process_id id() const { return SURFACE_ACOUSTIC_PHONON_SCATTERING; }

    protected:

      /** @brief Returns all possible final scattering states for a carrier with initial kinetic energy 'kin_energy' for an edge with vertices v1 and v2
       *
       * @param device     A ViennaSHE device object describing the physical device
       * @param controller The simulator controller
       * @param v1         First vertex of the edge
       * @param v2         Second vertex of the edge
       * @param kinetic_energy  Kinetic energy of the particle before scattering
       * @param ctype           Carrier type, either electrons or holes
       *
       * @return A vector describing all possible states. For ionized_impurity_scattering we use the approximation 'incoming energy equal final energy'.
       */
      template <typename ElementType>
      scatter_processes_type get(DeviceType const & device,
                                 ControllerType const & controller,
                                 ElementType const & elem,
                                 double kinetic_energy,
                                 viennashe::carrier_type_id ctype) const
      {
        scatter_processes_type result(1);

        // Final energy is always equal to initial energy:
        result[0].initial_energy(kinetic_energy);
        result[0].final_energy(kinetic_energy);
        result[0].rate(getScatteringRate(device, controller, elem, kinetic_energy, ctype));

        return result;
      }

      template <typename DeviceType, typename ControllerType, typename ElementType>
      double getScatteringRate(DeviceType const & device,
          ControllerType const & controller,
          ElementType const & elem,
          double kinetic_energy,
          viennashe::carrier_type_id ctype) const
      {
        const double temp = device.get_lattice_temperature(elem);
        const double Ep   = this->get_electric_field_n(device, elem);

        double rate = this->getScatteringRate(temp, Ep, ctype);

        // check if the rate is numerically zero
        if ( std::abs(rate) < std::numeric_limits<double>::epsilon() ) rate = 0.0; // assign zero to ensure matrix sparsity

        return rate;
      }

      template < typename DeviceType >
      double get_electric_field_n(DeviceType const & device,
          typename viennagrid::result_of::ncell<typename DeviceType::mesh_type::config_type, 0 > ::type const & vt) const
      {
        typedef typename DeviceType::mesh_type MeshType;
        typedef typename MeshType::config_type Config;
        typedef typename Config::cell_tag CellTag;
        typedef typename viennagrid::result_of::point<Config>::type PointType;

        PointType E;
        typename ElectricFieldAccessor::value_type Eacc = _Efield(vt);
        for ( std::size_t i = 0; i < PointType::dim; i++) E[i] = Eacc[i];

        PointType n = device.vector_to_next_insulator(vt);
        const double distance = viennagrid::norm_2(n);
        if ( distance > 0) n /= distance; // normalise

        return viennagrid::inner_prod(E, n); // projection on the vector towards the interface
      }

      template < typename DeviceType >
      double get_electric_field_n(DeviceType const & device,
          typename viennagrid::result_of::ncell<typename DeviceType::mesh_type::config_type, 1 > ::type const & edge) const
      {
        typedef typename DeviceType::mesh_type MeshType;
        typedef typename MeshType::config_type Config;
        typedef typename Config::cell_tag CellTag;
        typedef typename viennagrid::result_of::point<Config>::type PointType;
        typedef typename viennagrid::result_of::ncell<typename DeviceType::mesh_type::config_type, 0 > ::type VertexType;

        const double Emag = _Efield(edge); // field along the edge
        VertexType const & vt1 = viennagrid::ncells < 0 > (edge)[0];
        VertexType const & vt2 = viennagrid::ncells < 0 > (edge)[1];

        PointType n = 0.5 * ( device.vector_to_next_insulator(vt1) + device.vector_to_next_insulator(vt2) );
        const double distance = viennagrid::norm_2(n);
        if ( distance > 0) n /= distance; // normalise

        PointType E = (vt2.point() - vt1.point()) ;
        if (viennagrid::norm_2(E) > 0) E /= viennagrid::norm_2(E);
        E *= Emag;
        return viennagrid::inner_prod(E, n); // projection on the vector towards the interface
      }

      double getSurfaceInversionLayerWidth(const double T, const double Ep, viennashe::carrier_type_id ctype) const
      {
        const double q     = viennashe::physics::constants::q;
        const double kB    = viennashe::physics::constants::kB;
        const double delta = _params.get_potential_well_strength(ctype);
        const double hbar  = viennashe::physics::constants::hbar;
        const double mstar = viennashe::materials::si::dos_effective_mass(ctype);

        if (Ep == 0.0) return 0.0; // avoid division by zero ... in the limit of Ep -> 0 this is 0 anyway
        return q / (3.0 * kB * T / (2.0 * Ep) + delta * std::pow((hbar * q * hbar * q) / (mstar * Ep), 1.0 / 3.0)); // 1/m
      }

      double getScatteringRate(const double T, const double Ep, viennashe::carrier_type_id ctype) const
      {
        const double kB   = viennashe::physics::constants::kB;
        const double hbar = viennashe::physics::constants::hbar;
        const double ul   = _params.get_longitudinal_sound_velocity(ctype);
        const double E1   = _params.get_deformation_potential(ctype);
        const double rho  = _params.get_mass_density(ctype);
        const double pi   = viennashe::math::constants::pi;

        if (Ep == 0.0) return 0.0; // avoid division by zero

        // MUST NOT include density of states and/or N_op
        double scattering_rate = 2.0 * pi * kB * T * E1 * E1 / (hbar * ul * ul * rho);
        scattering_rate = scattering_rate * this->getSurfaceInversionLayerWidth(T, Ep, ctype);

        //log::debug<log_acoustic_phonon_scattering>() << "acoustical phonon: " << scattering_rate << std::endl;

        return _params.get_fit_factor() * scattering_rate;
      }

    private:

      viennashe::materials::surface_acoustic_phonon_scattering_parameters _params;
      ElectricFieldAccessor const & _Efield;


    };

  } //namespace she
} //namespace viennashe


#endif

#endif
