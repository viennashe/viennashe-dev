#ifndef VIENNASHE_SHE_POSTPROC_MACROSCOPIC_HPP
#define VIENNASHE_SHE_POSTPROC_MACROSCOPIC_HPP

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

/** @file viennashe/she/postproc/macroscopic.hpp
    @brief Writer for arbitrary macroscopic quantities
*/


namespace viennashe
{
  /** @brief Writes the provided macroscopic quantity to the container provided
   *
   * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param quantity       An accessor for the macroscopic quantity
   * @param quantity_container        The container to be filled with values
   * @tparam ContainerType   A container type, for example: std::vector or std::deque
   */
  template <typename DeviceType,
            typename MacroscopicQuantityAccessor,
            typename ContainerType>
  void write_macroscopic_quantity_to_container(DeviceType const & device,
                                               MacroscopicQuantityAccessor const & quantity,
                                               ContainerType & quantity_container)
  {
    typedef typename DeviceType::mesh_type              MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type    CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type       CellIterator;

    MeshType const & mesh = device.mesh();

    CellContainer cells(mesh);
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {
      quantity_container.at(static_cast<std::size_t>(cit->id().get())) = quantity(*cit);
    }
  }

} //namespace viennashe

#endif
