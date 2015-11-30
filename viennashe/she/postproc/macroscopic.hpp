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
#include <vector>

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
  template <typename DeviceType, typename MacroscopicQuantityAccessor>
  void write_macroscopic_quantity_to_quantity_field(DeviceType const & device,
                                                    MacroscopicQuantityAccessor const & quantity,
                                                    viennagrid_quantity_field field)
  {
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      std::vector<viennagrid_numeric> value = quantity(*cit);
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_quantity_field_value_set(field, *cit, &(value[0])));
    }
  }

} //namespace viennashe

#endif
