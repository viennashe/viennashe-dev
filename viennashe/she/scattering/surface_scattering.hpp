#ifndef VIENNASHE_SHE_SCATTERING_SURFACE_SCATTERING_HPP
#define VIENNASHE_SHE_SCATTERING_SURFACE_SCATTERING_HPP

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

#include "viennashe/accessors.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/scattering/common.hpp"

// viennagrid
#include "viennagrid/forwards.hpp"
#include "viennagrid/algorithm/inner_prod.hpp"
#include "viennagrid/algorithm/norm.hpp"

/** @file viennashe/she/scattering/surface_scattering.hpp
    @brief Implements the surface scattering processes using a phenomenological description (Lombardi).
 */

namespace viennashe
{
  namespace she
  {

    /** @brief Surface scattering process
     *
     *  Calculates surface scattering based on an estimate of the wave-function of the carriers using the electric field
     */
    template <typename DeviceType, typename ElectricFieldAccessor>
    class surface_scattering : public scattering_base<DeviceType>
    {
        typedef scattering_base<DeviceType>                     base_type;
        typedef typename base_type::FacetType                   FacetType;
        typedef typename base_type::CellType                    CellType;

      public:
        typedef typename base_type::scatter_processes_type      scatter_processes_type;
        typedef scatter_processes_type                          value_type;

      explicit surface_scattering(DeviceType const & device,
                                  viennashe::config const & conf,
                                  ElectricFieldAccessor const & Efield)
       : base_type(device, conf), params_(conf.scattering().surface()), _Efield(Efield) { }

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

        scatter_process_id id() const { return SURFACE_SCATTERING; }

    protected:

      /** @brief Returns all possible final scattering states for a carrier with initial kinetic energy 'kin_energy' for an edge with vertices v1 and v2
       *
       * @param elem       The element for which to obtain the scattering process
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
        result[0].rate(getScatteringRate(elem, kinetic_energy, ctype));

        return result;
      }


      template <typename ElementType>
      double getScatteringRate(ElementType const & elem,
                               double kinetic_energy,
                               viennashe::carrier_type_id ctype) const
      {
        (void)kinetic_energy; //avoid unused parameter warnings
        const double temp = base_type::device_.get_lattice_temperature(elem);
        const double Ep   = this->get_electric_field_n(elem);

        const double totaldoping = base_type::device_.get_doping_n(elem) + base_type::device_.get_doping_p(elem);

        double rate = this->getScatteringRate(temp, totaldoping, Ep, ctype);

        // check if the rate is numerically zero
        if ( std::abs(rate) < std::numeric_limits<double>::epsilon() ) rate = 0.0; // assign zero to ensure matrix sparsity

        return rate;
      }

      double get_electric_field_n(CellType const &) const
      {
        throw std::runtime_error("get_charged_trap_density(): TODO");
        /*
        typedef typename DeviceType::mesh_type                         MeshType;
        typedef typename viennagrid::result_of::point<MeshType>::type  PointType;

        PointType n = base_type::device_.vector_to_next_insulator(vt);
        const double distance = viennagrid::norm_2(n);
        if ( distance > 0) n /= distance; // normalise

        if ( distance > params_.cutoff_distance() ) return 0.0;

        PointType E;
        typename ElectricFieldAccessor::value_type Eacc = _Efield(vt);
        for ( std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); i++) E[i] = Eacc[i];

        return viennagrid::inner_prod(E, n); // projection on the vector towards the interface
        */
      }

      double get_electric_field_n(FacetType const & /*edge*/) const
      {
        throw std::runtime_error("get_charged_trap_density(): TODO");
        /*
        typedef typename DeviceType::mesh_type                         MeshType;
        typedef typename viennagrid::result_of::point<MeshType>::type  PointType;

        VertexType const & vt1 = viennagrid::vertices(edge)[0];
        VertexType const & vt2 = viennagrid::vertices(edge)[1];

        PointType n = 0.5 * ( base_type::device_.vector_to_next_insulator(vt1) + base_type::device_.vector_to_next_insulator(vt2) );
        const double distance = viennagrid::norm_2(n);
        if ( distance > 0) n /= distance; // normalise

        if ( distance > params_.cutoff_distance() ) return 0.0;

        const double Emag = _Efield(edge); // field along the edge

        PointType E = (viennagrid::point(vt2) - viennagrid::point(vt1)) ;
        if (viennagrid::norm_2(E) > 0) E /= viennagrid::norm_2(E);
        E *= Emag;
        return viennagrid::inner_prod(E, n); // projection on the vector towards the interface
        */
      }

      double getScatteringRate(const double T, const double totaldoping, const double Ep, viennashe::carrier_type_id ctype) const
      {
        (void)T; //avoid unused parameter warnings
        //const double kB   = viennashe::physics::constants::kB;
        //const double hbar = viennashe::physics::constants::hbar;
        const double pi   = viennashe::math::constants::pi;

        const double dos_mass = viennashe::materials::si::dos_effective_mass(ctype);

        const double prefactor = ( std::sqrt(2.0) / std::pow((2.0*pi),3.0) ) * std::pow(dos_mass, 1.5) * 2.0/std::sqrt(pi);

        if (Ep <= 0.0) return 0.0; // avoid division by zero and Nan

        double rate =   1.0 / ( params_.first_factor(ctype) / Ep
                      + params_.second_factor(ctype) * std::pow(totaldoping, 0.125) * 3.0 / std::pow(Ep, 1.0/3.0) )
                      + params_.third_factor(ctype) * Ep * Ep;
        rate *= 4.0 * pi;

        return params_.get_fit_factor(ctype) * prefactor * rate;
      }

    private:

      surface_scattering_parameters params_;
      ElectricFieldAccessor const & _Efield;


    };

  } //namespace she
} //namespace viennashe



#endif
