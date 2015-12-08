#ifndef VIENNASHE_SHE_POSTPROC_EXPANSION_ORDER_HPP
#define VIENNASHE_SHE_POSTPROC_EXPANSION_ORDER_HPP

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
#include "viennashe/she/postproc/macroscopic.hpp"

/** @file viennashe/she/postproc/expansion_order.hpp
    @brief Provides an accessor for the expansion order in (x,H)-space
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Accessor class providing the average expansion order inside the device */
    template <typename SHEQuantity>
    class average_expansion_order_wrapper
    {
      public:

        typedef double value_type;

        average_expansion_order_wrapper(SHEQuantity const & quan,
                                        double energy_start,
                                        double energy_end) : quan_(quan),
                                                             energy_start_(energy_start),
                                                             energy_end_(energy_end) {}

        std::vector<value_type> operator()(viennagrid_element_id cell) const
        {
          double sum_order = 0;   //summation of all expansion orders
          long num_energies = 0;  //no of energy points considered

          for (std::size_t index_H=1; index_H<quan_.get_value_H_size() - 1; ++index_H)
          {
            double kinetic_energy = quan_.get_kinetic_energy(cell, index_H);

            if (kinetic_energy >= energy_start_ && kinetic_energy <= energy_end_)
            {
              sum_order += quan_.get_expansion_order(cell, index_H);
              ++num_energies;
            }
          }

          std::vector<value_type> ret(3);
          ret[0] = sum_order / num_energies;
          return ret;
        }

      private:
        SHEQuantity quan_;
        double energy_start_;
        double energy_end_;
    };

    /** @brief Convenience function for writing the average expansion order to the container provided.
     *
     * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan           The SHE quantities. Passed for compatibility with other write_XYZ_to_domain() routines.
     * @param container      The container to write to
     * @param energy_start   (Optional) Lower bound for the kinetic energy (in eV) range to be considered for the averaging the SHE order
     * @param energy_end     (Optional) Upper bound for the kinetic energy (in eV) range to be considered for the averaging the SHE order
     */
    template <typename SHEQuantity>
    void write_average_expansion_order_to_container(viennashe::device const & device,
                                                    SHEQuantity const & quan,
                                                    viennagrid_quantity_field field,
                                                    double energy_start = 0.0,
                                                    double energy_end = 1.0)
    {
      average_expansion_order_wrapper<SHEQuantity> wrapper(quan, energy_start, energy_end);

      viennashe::write_macroscopic_quantity_to_quantity_field(device, wrapper, field);
    }

  } //namespace she
} //namespace viennashe

#endif
