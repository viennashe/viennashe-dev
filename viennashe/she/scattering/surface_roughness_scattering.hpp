#ifndef VIENNASHE_SHE_SCATTERING_SURFACE_ROUGHNESS_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_SURFACE_ROUGHNESS_SCATTERING_HPP

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
#include "viennagrid/forwards.hpp"
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
    class surface_roughness_scattering : public scattering_base<DeviceType, ControllerType>
    {
        typedef scattering_base<DeviceType, ControllerType>     base_type;
        typedef typename base_type::VertexType                  VertexType;
        typedef typename base_type::EdgeType                    EdgeType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

      explicit surface_roughness_scattering(const viennashe::materials::surface_roughness_scattering_parameters & params,
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

        scatter_process_id id() const { return SURFACE_ROUGHNESS_SCATTERING; }

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
        const double Ep   = this->get_electric_field_n(device, elem, ctype);

        //std::cout << "  the rate is   " << this->getScatteringRate(temp, Ep, ctype) << std::endl << std::endl;

        return this->getScatteringRate(temp, Ep, ctype);
      }

      template < typename DeviceType>
      double get_electric_field_n(DeviceType const & device,
          typename viennagrid::result_of::ncell<typename DeviceType::mesh_type::config_type, 0 > ::type const & vt,
          viennashe::carrier_type_id ctype) const
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

        if ( this->_params.get_correlation_length(tag) < distance )
        {
          if ( distance > 0) n /= distance; // normalise
          //std::cout << vt.id() << " => " << distance << " with (" << n << ")  |  " << E << " projected " << viennagrid::inner_prod(E, n) << std::endl;
          //std::cout << " SURFACE ROUGHNESS SCATTERING Ep = " << viennagrid::inner_prod(E, n) << " AT " << vt << std::endl;
          return viennagrid::inner_prod(E, n); // projection on the vector towards the interface
        }

        return 0.0; // outside of the correlation length the scattering rate is zero
      }

      template < typename DeviceType>
      double get_electric_field_n(DeviceType const & device,
          typename viennagrid::result_of::ncell<typename DeviceType::mesh_type::config_type, 1 > ::type const & edge,
          viennashe::carrier_type_id ctype) const
      {
        // TODO: CHECK THIS ... is this a good idea?
        return 0.5 * (this->get_electric_field_n(device, viennagrid::ncells < 0 > (edge)[0], tag) +
                      this->get_electric_field_n(device, viennagrid::ncells < 0 > (edge)[1], tag));
      }

      double getEffectiveField(const double T, const double Ep, viennashe::carrier_type_id ctype) const
      {
        return Ep; // TODO: maybe introduce a fit formula here
      }

      double getScatteringRate(const double T, const double Ep, viennashe::carrier_type_id ctype) const
      {
        const double hbar = viennashe::physics::constants::hbar;
        const double pi   = viennashe::math::constants::pi;
        const double q    = viennashe::physics::constants::q;
        const double drms = _params.get_rms_height(ctype);
        const double L    = _params.get_correlation_length(ctype);
        const double esi  = viennashe::materials::si::permittivity();
        const double Eeff = this->getEffectiveField(T, Ep, ctype);

        //double scattering_rate = 2.0 * pi * q * q / hbar * (esi / (esi+eox)) * (esi / (esi+eox)) * drms * drms * L * L; // TODO: Use improved rate ?
        double scattering_rate = 2.0 * pi * q * q / hbar * (esi * esi) * drms * drms * L * L;
        scattering_rate = scattering_rate * Eeff * Eeff;

        return _params.get_fit_factor() * scattering_rate;
      }

    private:

      viennashe::materials::surface_roughness_scattering_parameters _params;
      ElectricFieldAccessor const & _Efield;

    };

  } //namespace she
} //namespace viennashe

#endif

#endif /* VIENNASHE_SHE_SCATTERING_SURFACE_ROUGHNESS_SCATTERING_HPP */

