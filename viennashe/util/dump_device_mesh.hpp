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

#include "viennagrid/forwards.hpp"

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
      typedef typename viennagrid::result_of::point<MeshType>::type      PointType;
      typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;
      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type        CellContainer;

      typedef typename viennagrid::result_of::const_vertex_range<CellType>::type     VertexOnCellContainer;
      typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type  VertexOnCellIterator;

      if (vertices == 0) throw std::invalid_argument("vertices = NULL !");
      if (cells == 0)    throw std::invalid_argument("cells = NULL !");

      VertexContainer grid_vertices(device.mesh());
      CellContainer   grid_cells(device.mesh());
      const size_t    dim = PointType::dim;

      // Dump sizes
      num_vertices = static_cast<IndexT>(grid_vertices.size());
      num_cells    = static_cast<IndexT>(grid_cells.size());

      // Dump vertices (iterate the C-way to have the indices)
      for (std::size_t i = 0; i < grid_vertices.size(); ++i)
      {
        for (std::size_t j = 0; j < dim; ++j) vertices[i][j] = viennagrid::point(grid_vertices[i])[j];
      }

      // Dump cells (iterate the C-way to have the indices)
      for (std::size_t i = 0; i < grid_cells.size(); ++i)
      {
        VertexOnCellContainer vertices_on_cell(grid_cells[i]);
        std::size_t j = 0;
        for (VertexOnCellIterator vit = vertices_on_cell.begin(); vit != vertices_on_cell.end(); ++vit, ++j)
        {
          cells[i][j] = static_cast<IndexT>(vit->id().get());
        }
      }

    } // dump_mesh

  }
} // namespace viennashe


#endif	/* VIENNASHE_UTIL_DUMP_DEVICE_MESH_HPP */

