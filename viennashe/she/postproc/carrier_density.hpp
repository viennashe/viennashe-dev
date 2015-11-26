#ifndef VIENNASHE_SHE_CARRIER_DENSITY_HPP
#define VIENNASHE_SHE_CARRIER_DENSITY_HPP

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


#include <math.h>
#include <fstream>
#include <iostream>

#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"
#include "viennashe/config.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

/** @file viennashe/she/postproc/carrier_density.hpp
    @brief Provides an accessor for the carrier density
*/

namespace viennashe
{
  namespace she
  {

    namespace detail
    {
      /** @brief An accessor for the carrier density in the device by reference  */
      template <typename SHEQuantity>
      class carrier_density_wrapper_by_reference
      {
          typedef viennashe::config::dispersion_relation_type      dispersion_relation_type;

        public:

          typedef double value_type;

          carrier_density_wrapper_by_reference(viennashe::config const & conf,
                                  SHEQuantity const & quan,
                                  double energy_start = 0.0,
                                  double energy_end = 1.0) : conf_(conf), quan_(quan), energy_start_(energy_start), energy_end_(energy_end) {}

          carrier_density_wrapper_by_reference(carrier_density_wrapper_by_reference const & o)
            : quan_(o.quan_), conf_(o.conf_), energy_start_(o.energy_start_), energy_end_(o.energy_end_)
          { }

          carrier_density_wrapper_by_reference(carrier_density_wrapper_by_reference const & o, SHEQuantity const & quan)
            : conf_(o.conf_), quan_(quan), energy_start_(o.energy_start_), energy_end_(o.energy_end_)
          { }

          template <typename ElementType>
          value_type operator()(ElementType const & elem) const
          {
            value_type density = 0;

            typename viennashe::config::dispersion_relation_type dispersion = conf_.dispersion_relation(quan_.get_carrier_type_id());

            if (quan_.get_unknown_mask(elem))
            {
              for (std::size_t index_H=1; index_H < quan_.get_value_H_size() - 1; ++index_H)
              {
                const long unknown_index = quan_.get_unknown_index(elem, index_H);
                if (unknown_index < 0)
                  continue;

                const double energy_lower = std::max( (quan_.get_kinetic_energy(elem, index_H - 1) + quan_.get_kinetic_energy(elem, index_H)) / 2.0,
                                                energy_start_);
                const double energy_upper = std::min( (quan_.get_kinetic_energy(elem, index_H + 1) + quan_.get_kinetic_energy(elem, index_H)) / 2.0,
                                                  energy_end_);

                if (energy_upper >= 0 && energy_lower < energy_upper)
                {
                  double height = box_height(quan_, elem, index_H);
                  switch (conf_.she_discretization_type())
                  {
                  case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                    density += quan_.get_values(elem, index_H)[0] * averaged_density_of_states(quan_, dispersion, elem, index_H) * height;
                    break;
                  case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                    density += quan_.get_values(elem, index_H)[0] * height;
                    break;
                  default: throw std::runtime_error("carrier_density_wrapper_by_reference::operator(): Unknown SHE discretization type!");
                  }
                }
              } // for index_H
            }
            else // If there's a boundary condition
            {
              for (std::size_t index_H=1; index_H < quan_.get_value_H_size() - 1; ++index_H)
              {
                const double energy_lower = std::max( (quan_.get_kinetic_energy(elem, index_H - 1) + quan_.get_kinetic_energy(elem, index_H)) / 2.0,
                                                energy_start_);
                const double energy_upper = std::min( (quan_.get_kinetic_energy(elem, index_H + 1) + quan_.get_kinetic_energy(elem, index_H)) / 2.0,
                                                  energy_end_);
                if (energy_upper >= 0 && energy_lower < energy_upper)
                {
                  double height = box_height(quan_, elem, index_H);
                  switch (conf_.she_discretization_type())
                  {
                  case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                    density += quan_.get_boundary_value(elem, index_H) * averaged_density_of_states(quan_, dispersion, elem, index_H) * height;
                    break;
                  case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                    density += quan_.get_boundary_value(elem, index_H) * height;
                    break;
                  default: throw std::runtime_error("carrier_density_wrapper_by_reference::operator(): Unknown SHE discretization type!");
                  }
                }
              } // for index_H
            }

            if ( viennashe::util::is_Inf(density) ) viennashe::log::warn() << "* carrier_density_wrapper()::operator(): WARNING: density is inf " << elem << std::endl;
            if ( viennashe::util::is_NaN(density) ) viennashe::log::warn() << "* carrier_density_wrapper()::operator(): WARNING: density is nan " << elem << std::endl;

            return density * dispersion.symmetry_factor();
          }

        private:
          viennashe::config conf_; // NO REFERENCE!
          SHEQuantity const & quan_;
          double energy_start_;
          double energy_end_;
      };
    } // namespace detail

    /** @brief An accessor for the carrier density in the device */
    template <typename SHEQuantity>
    class carrier_density_wrapper
    {
        typedef viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef viennashe::she::detail::carrier_density_wrapper_by_reference<SHEQuantity> carrier_density_by_ref_type;

      public:

        typedef typename carrier_density_by_ref_type::value_type value_type;

        carrier_density_wrapper(viennashe::config const & conf,
                                SHEQuantity const & quan,
                                double energy_start = 0.0,
                                double energy_end = 1.0)
          : quan_(quan /* COPY ALL VALUES */),
            density_impl_by_ref_(conf, quan_ /* REFERENCE COPIED VALUES */, energy_start, energy_end) { }

        carrier_density_wrapper(carrier_density_wrapper const & o)
          : quan_(o.quan_ /* COPY ALL VALUES */),
            density_impl_by_ref_(o.density_impl_by_ref_, quan_ /* REFERENCE COPIED VALUES */)
        { }

        template <typename ElementType>
        value_type operator()(ElementType const & elem) const
        {
          return this->density_impl_by_ref_(elem);
        }

      private:
        SHEQuantity quan_;
        carrier_density_by_ref_type density_impl_by_ref_;
    };


    /** @brief Computes the carrier density in the device and writes the result to a vector. Used during Gummel iteration. */
    template <typename DeviceType,
              typename SHEQuantity,
              typename VectorType>
    void compute_carrier_density_vector(DeviceType const & device,
                                        SHEQuantity const & quan,
                                        viennashe::config::dispersion_relation_type const & dispersion,
                                        VectorType & carrier)
    {
      typedef typename DeviceType::mesh_type              MeshType;
      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type    CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type       CellIterator;

      MeshType const & mesh = device.mesh();

      carrier_density_wrapper<SHEQuantity> carrier_dens(quan, dispersion);

      //carrier concentration is computed for even unknowns, which are associated with vertices:
      CellContainer cells(mesh);
      for (CellIterator cit = cells.begin();
            cit != cells.end();
            ++cit)
      {
        long carrier_index = quan.get_unknown_index(*cit);

        //only compute in the interior
        if (carrier_index >= 0)
          carrier[carrier_index] = carrier_dens(*cit);
      }
    }


    /** @brief Convenience function for writing the average expansion order to the container provided.
     *
     * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan           The SHE quantities. Passed for compatibility with other write_XYZ_to_domain() routines.
     * @param dispersion     The dispersion relation
     * @param container      Container for the quantity
     * @param energy_start   (Optional) Lower bound for the kinetic energy (in eV) range to be considered for the averaging the SHE order
     * @param energy_end     (Optional) Upper bound for the kinetic energy (in eV) range to be considered for the averaging the SHE order
     */
    template <typename DeviceType,
              typename SHEQuantity,
              typename ContainerType>
    void write_carrier_density_to_container(DeviceType const & device,
                                            SHEQuantity const & quan,
                                            viennashe::config::dispersion_relation_type const & dispersion,
                                            ContainerType & container,
                                            double energy_start = 0.0,
                                            double energy_end = 1.0)
    {
      carrier_density_wrapper<SHEQuantity> wrapper(quan, dispersion, energy_start, energy_end);

      viennashe::write_macroscopic_quantity_to_container(device, wrapper, container);
    }



  } //namespace she
} //namespace viennashe

#endif
