#ifndef VIENNASHE_SHE_CARRIER_VELOCITY_HPP
#define VIENNASHE_SHE_CARRIER_VELOCITY_HPP

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
#include "viennashe/physics/constants.hpp"
#include "viennashe/math/constants.hpp"
#include "viennashe/math/integrator.hpp"
#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/she/assemble_common.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/she/postproc/current_density.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"


/** @file viennashe/she/postproc/carrier_velocity.hpp
    @brief Provides an accessor for computing the average carrier drift velocity in each (semiconductor) cell.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Accessor class providing the carrier velocity inside the device */
    template <typename DeviceType,
              typename SHEQuantity>
    class carrier_velocity_wrapper
    {
      public:
        typedef std::vector<double>       value_type;

        carrier_velocity_wrapper(DeviceType        const & device,
                                 viennashe::config const & conf,
                                 SHEQuantity       const & quan)
          : device_(device), quan_(quan),
            current_density_(device, conf, quan), carrier_density_(conf, quan) {}

        /** @brief Functor interface returning the average carrier drift velocity at the provided cell */
        template <typename CellT>
        value_type operator()(CellT const & cell) const
        {
          typedef typename viennagrid::result_of::point<CellT>::type          PointType;

          std::vector<double> carrier_velocity(3);

          if (!viennashe::materials::is_semiconductor(device_.get_material(cell)))
            return carrier_velocity;

          double value_Y_00 = viennashe::math::SphericalHarmonic(0,0)(0.0, 0.0);
          double polarity = (quan_.get_carrier_type_id() == ELECTRON_TYPE_ID) ? -1.0 : 1.0;

          carrier_velocity = current_density_(cell);
          double carrier_density = carrier_density_(cell);

          for (std::size_t i=0; i<static_cast<std::size_t>(PointType::dim); ++i)
            carrier_velocity[i] = -polarity * value_Y_00 * carrier_velocity[i] / carrier_density / viennashe::physics::constants::q;  //cf. current_density_wrapper on these additional prefactors

          return carrier_velocity;
        }

      private:

        DeviceType  const & device_;
        SHEQuantity const & quan_;

        current_density_wrapper<DeviceType, SHEQuantity> current_density_;
        carrier_density_wrapper<SHEQuantity>             carrier_density_;
    };


    /** @brief Convenience function for writing the average expansion order to the container provided.
     *
     * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param conf             The simulator configuration
     * @param quan             The SHE quantity (electron or hole distribution function) used for the calculation
     * @param container        The container the values should be written to
    */
    template <typename DeviceType,
              typename SHEQuantity,
              typename ContainerType>
    void write_carrier_velocity_to_container(DeviceType const & device,
                                             viennashe::config const & conf,
                                             SHEQuantity const & quan,
                                             ContainerType & container)
    {
      carrier_velocity_wrapper<DeviceType, SHEQuantity> velocity_wrapper(device, conf, quan);

      viennashe::write_macroscopic_quantity_to_container(device, velocity_wrapper, container);
    }


  } //namespace she
} //namespace viennashe

#endif
