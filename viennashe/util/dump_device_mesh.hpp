#ifndef VIENNASHE_UTIL_DUMP_DEVICE_MESH_HPP
#define	VIENNASHE_UTIL_DUMP_DEVICE_MESH_HPP
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

#include "viennashe/forwards.h"
#include "viennashe/device.hpp"

#include "viennagrid/viennagrid.h"

namespace viennashe
{
  namespace util
  {

    template < typename DeviceT, typename IndexT >
    void dump_mesh(DeviceT const & device, double ** vertices, IndexT & num_vertices,
                   IndexT ** cells, IndexT & num_cells
                  )
    {
      typedef typename DeviceT::mesh_type MeshType;

      if (vertices == 0) throw std::invalid_argument("vertices = NULL !");
      if (cells == 0)    throw std::invalid_argument("cells = NULL !");

      viennagrid_element_id *vertices_begin, *vertices_end;
      viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

      viennagrid_dimension cell_dim;
      viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

      viennagrid_element_id *cells_begin, *cells_end;
      viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

      // Dump sizes
      num_vertices = IndexT(vertices_end - vertices_begin);
      num_cells    = IndexT(   cells_end -    cells_begin);

      viennagrid_dimension geo_dim;
      viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim);

      // Dump vertices (iterate the C-way to have the indices)
      for (std::size_t i = 0; i < num_vertices; ++i)
      {
        viennagrid_numeric *coords;
        viennagrid_mesh_vertex_coords_get(device.mesh(), vertices_begin[i], &coords);

        for (std::size_t j = 0; j < std::size_t(geo_dim); ++j)
          vertices[i][j] = coords[j];
      }

      // Dump cells (iterate the C-way to have the indices)
      for (std::size_t i = 0; i < num_cells; ++i)
      {
        viennagrid_element_id *vertices_on_cell_begin, *vertices_on_cell_end;
        viennagrid_element_boundary_elements(device.mesh(), cells_begin[i], 0, &vertices_on_cell_begin, &vertices_on_cell_end);

        std::size_t j = 0;
        for (viennagrid_element_id *vit = vertices_on_cell_begin; vit != vertices_on_cell_end; ++vit, ++j)
        {
          cells[i][j] = IndexT(viennagrid_index_from_element_id(*vit));
        }
      }

    } // dump_mesh

  }
} // namespace viennashe


#endif	/* VIENNASHE_UTIL_DUMP_DEVICE_MESH_HPP */

