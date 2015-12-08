#ifndef VIENNASHE_SHE_POSTPROC_POTENTIAL_HPP
#define VIENNASHE_SHE_POSTPROC_POTENTIAL_HPP

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

#include "viennagrid/mesh/mesh.hpp"

#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"

/** @file viennashe/postproc/potential.hpp
    @brief Writer for the potential
*/

namespace viennashe
{
  /** @brief Writes the provided potential to the container provided
   *
   * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param pot            An accessor for the potential
   * @param container        The container to be filled with values
   * @tparam ContainerType   A container type, for example: std::vector or std::deque
   */
  template <typename PotentialAccessor,
            typename ContainerType>
  void write_potential_to_quantity_field(viennashe::device const & device,
                                         PotentialAccessor const & pot,
                                         viennagrid_quantity_field field)
  {
    viennashe::write_macroscopic_quantity_to_quantity_field(device, pot, field);
  }

} //namespace viennashe

#endif
