#ifndef VIENNASHE_SHE_CARRIER_KINETIC_ENERGY_HPP
#define VIENNASHE_SHE_CARRIER_KINETIC_ENERGY_HPP

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



#include <math.h>
#include <fstream>
#include <iostream>

#include "viennagrid/mesh/mesh.hpp"

#include "viennashe/forwards.h"
#include "viennashe/config.hpp"
#include "viennashe/physics/constants.hpp"

#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"

/** @file viennashe/she/postproc/carrier_energy.hpp
    @brief Provides an accessor for the average carrier energy
*/

namespace viennashe
{
  namespace she
  {
    namespace detail
    {
      class energy_integrator_for_box
      {
          typedef viennashe::config::dispersion_relation_type      dispersion_relation_type;
        public:
          energy_integrator_for_box(dispersion_relation_type const & dispersion,
                                    double f_00_value) : dispersion_(dispersion), f_00_(f_00_value) {}

          double operator()(double kinetic_energy) const
          {
            return  kinetic_energy
                    * dispersion_.symmetry_factor()     //in order to get consistent results when dividing by the carrier density
                    * f_00_;
          }

        private:
          dispersion_relation_type dispersion_;
          double f_00_;
      };
    } // namespace detail


    /** @brief An accessor for the average carrier energy at each point inside the device */
    template <typename SHEQuantity>
    class carrier_energy_wrapper
    {
        typedef viennashe::config::dispersion_relation_type      dispersion_relation_type;

    public:

      typedef double value_type;

      carrier_energy_wrapper(viennashe::config const & conf,
                             SHEQuantity const & quan)
          : conf_(conf), quan_(quan) {}

      template <typename CellType>
      value_type operator()(CellType const & cell) const
      {
        typedef typename detail::carrier_density_wrapper_by_reference<SHEQuantity>     carrier_density_type;
        carrier_density_type carrier_density_(conf_, quan_);

        double avg_energy = 0;
        const double density = carrier_density_(cell);

        typename viennashe::config::dispersion_relation_type dispersion = conf_.dispersion_relation(quan_.get_carrier_type_id());

        if (density <= 0)
          return 0;

        if (quan_.get_unknown_mask(cell))
        {
          for (std::size_t index_H=1; index_H < quan_.get_value_H_size() - 1; ++index_H)
          {
            if ( quan_.get_unknown_index(cell, index_H) < 0 )
              continue;

            const double energy_mid   = std::max(quan_.get_kinetic_energy(cell, index_H), 0.0);
            const double energy_lower = std::max( (quan_.get_kinetic_energy(cell, index_H - 1) + quan_.get_kinetic_energy(cell, index_H)) / 2.0,
                                                   0.0);
            const double energy_upper = (quan_.get_kinetic_energy(cell, index_H + 1) + quan_.get_kinetic_energy(cell, index_H)) / 2.0;

            if (energy_upper >= 0 && energy_lower < energy_upper)
            {
              double height = box_height(quan_, cell, index_H);
              switch (conf_.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                avg_energy += quan_.get_values(cell, index_H)[0] * energy_mid * averaged_density_of_states(quan_, dispersion, cell, index_H) * height;
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                avg_energy += quan_.get_values(cell, index_H)[0] * energy_mid * height;
                break;
              default: throw std::runtime_error("carrier_energy_wrapper::operator(): Unknown SHE discretization type!");
              }
            }
          } // index_H
        }
        else // Get the average carrier energy at the boundary (we've a boundary condition here)
        {
          for (std::size_t index_H=1; index_H < quan_.get_value_H_size() - 1; ++index_H)
          {
            const double energy_mid   = std::max(quan_.get_kinetic_energy(cell, index_H), 0.0);
            const double energy_lower = std::max( (quan_.get_kinetic_energy(cell, index_H - 1) + quan_.get_kinetic_energy(cell, index_H)) / 2.0, 0.0);
            const double energy_upper =(quan_.get_kinetic_energy(cell, index_H + 1) + quan_.get_kinetic_energy(cell, index_H)) / 2.0;

            if (energy_upper >= 0 && energy_lower < energy_upper)
            {
              double height = box_height(quan_, cell, index_H);
              switch (conf_.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                avg_energy += quan_.get_boundary_value(cell, index_H) * energy_mid * averaged_density_of_states(quan_, dispersion, cell, index_H) * height;
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                avg_energy += quan_.get_boundary_value(cell, index_H) * energy_mid * height;
                break;
              default: throw std::runtime_error("carrier_energy_wrapper::operator(): Unknown SHE discretization type!");
              }
            }
          } // index_H
        }

        return avg_energy / density;
      }

    private:
      viennashe::config conf_;
      SHEQuantity  quan_;
    };


    /** @brief Convenience function for writing the average kinetic carrier energy to the container provided.
     *
     * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan             The SHE quantities. Passed for compatibility with other write_XYZ_to_domain() routines.
     * @param conf             The simulator configuration
     * @param container        The container for storing the average carrier energy values
     */
    template <typename DeviceType,
              typename SHEQuantity,
              typename ContainerType>
    void write_kinetic_carrier_energy_to_container(DeviceType const & device,
                                                   viennashe::config const & conf,
                                                   SHEQuantity const & quan,
                                                   ContainerType & container)
    {
      viennashe::she::carrier_energy_wrapper<SHEQuantity> wrapper(conf, quan);

      viennashe::write_macroscopic_quantity_to_container(device, wrapper, container);
    }


  } //namespace she
} //namespace viennashe

#endif
