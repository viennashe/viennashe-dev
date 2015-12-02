#ifndef VIENNASHE_UTIL_TRANSLATE_QUANTITY_HPP
#define VIENNASHE_UTIL_TRANSLATE_QUANTITY_HPP

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

/** @file viennashe/util/dual_box_flux.hpp
    @brief Helper routines for projecting normal components of a vector-valued quantity defined on edges to vertices (e.g. for visualization purposes)
 */

// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/util/misc.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/math/linalg_util.hpp"

#include "viennashe/solvers/forwards.h"

namespace viennashe
{
  namespace util
  {

    inline std::vector<viennagrid_numeric>
    outer_cell_normal_at_facet(viennagrid_mesh mesh, viennagrid_element_id cell, viennagrid_element_id facet)
    {
      std::vector<viennagrid_numeric> centroid_cell(3);
      std::vector<viennagrid_numeric> centroid_facet(3);

      viennagrid_dimension cell_dim = viennagrid_topological_dimension_from_element_id(cell);

      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh,  cell, &(centroid_cell[0])));
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, facet, &(centroid_facet[0])));

      std::vector<viennagrid_numeric> ret(3);

      if (cell_dim == 1)
      {
        if (centroid_cell[0] < centroid_facet[0])
          ret[0] = 1.0;
        else
          ret[0] = -1.0;

        return ret;
      }


      viennagrid_element_id *facet_vertices_begin, *facet_vertices_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(mesh, facet, 0, &facet_vertices_begin, &facet_vertices_end));

      viennagrid_numeric *coords_0, *coords_1;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_vertex_coords_get(mesh, facet_vertices_begin[0], &coords_0));
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_vertex_coords_get(mesh, facet_vertices_begin[1], &coords_1));

      if (cell_dim == 2) // triangles and quadrilaterals
      {
        // create vector normal to facet: (x, y) -> (y, -x)
        ret[0] =    coords_1[1] - coords_0[1];
        ret[1] = - (coords_1[0] - coords_0[0]);
      }
      else if (cell_dim == 3) // tetrahedra and hexahedra
      {
        viennagrid_numeric *coords_2;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_vertex_coords_get(mesh, facet_vertices_begin[2], &coords_2));

        viennagrid_numeric direction_1[3] = { coords_1[0] - coords_0[0],
                                              coords_1[1] - coords_0[1],
                                              coords_1[2] - coords_0[2] };

        viennagrid_numeric direction_2[3] = { coords_2[0] - coords_0[0],
                                              coords_2[1] - coords_0[1],
                                              coords_2[2] - coords_0[2] };

        VIENNASHE_VIENNAGRID_CHECK(viennagrid_cross_prod(&(direction_1[0]), &(direction_2[0]), &ret[0]));
      }
      else
        throw std::runtime_error("outer_cell_normal_at_facet: invalid cell dimension!");

      // check direction: normal vector needs to point outwards, i.e. similar direction than connection cell_centroid to facet_centroid
      viennagrid_numeric inner_prod = 0;
      inner_prod += ret[0] * (centroid_facet[0] - centroid_cell[0]);
      inner_prod += ret[1] * (centroid_facet[1] - centroid_cell[1]);
      inner_prod += ret[2] * (centroid_facet[2] - centroid_cell[2]);

      viennagrid_numeric norm = std::sqrt(ret[0] * ret[0] + ret[1] * ret[1] + ret[2] * ret[2]);
      ret[0] /= (inner_prod < 0) ? -norm : norm;
      ret[1] /= (inner_prod < 0) ? -norm : norm;
      ret[2] /= (inner_prod < 0) ? -norm : norm;

      return ret;
    }


    /**
     * @brief Interpolates normal components of the flux defined on each facet to cells. Mostly used for visualization purposes.
     *
     * @param device           The device object (needed because it holds the Voronoi data)
     * @param cell             The cell to be interpolated to
     * @param facets           The set of facets to be considered (facets need to be the facets of the given cell!)
     * @param cell_setter      Functor for storing the interpolated flux vector
     * @param facet_access     Functor for accessing the normal components of the flux
     */
    template <typename DeviceT, typename CellSetterT, typename FacetAccessorT>
    void dual_box_flux_to_cell(DeviceT     const & device,
                               viennagrid_element_id cell,
                               CellSetterT       & cell_setter, FacetAccessorT  const & facet_access)
    {
      viennagrid_dimension geo_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim));

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *facets_begin, *facets_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(device.mesh(), cell, cell_dim - 1, &facets_begin, &facets_end));

      viennagrid_int facet_num = facets_end - facets_begin;

      viennashe::math::dense_matrix<double> M(geo_dim, geo_dim);
      std::vector<double>                   b(geo_dim);

      std::vector<std::vector<double> >  normals(facet_num);
      std::vector<double>     flux_contributions(facet_num);

      // for all edges connected with the current vertex
      std::size_t facet_ctr = 0;

      for (viennagrid_element_id *focit  = facets_begin;
                                  focit != facets_end;
                                ++focit, ++facet_ctr )
      {
        flux_contributions[facet_ctr] = facet_access(*focit);
        normals[facet_ctr]            = outer_cell_normal_at_facet(device.mesh(), cell, *focit); // vector along the edge pointing away from Vertex 'to'

        // flip orientation of flux contribution if global orientation is different:
        viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), *focit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));
        if (cells_on_facet_begin[0] != cell)
          flux_contributions[facet_ctr] *= -1.0;
      }

      // assemble mass matrix and rhs: M_i * V_i = rhs_i
      for ( std::size_t i = 0; i < static_cast<std::size_t>(geo_dim); ++i )
      {
        for ( std::size_t j = 0; j < static_cast<std::size_t>(geo_dim); ++j )
        {
          for ( std::size_t k = 0; k < normals.size(); ++k )
            M(i, j) += normals[k][i] * normals[k][j];
        }
        for ( std::size_t k = 0; k < normals.size(); ++k )
          b[i] += normals[k][i] * flux_contributions[k];
      }

      std::vector<double> to_value(geo_dim); // default: zero vector

      std::vector<double> result = viennashe::solvers::solve(M, b);
      for (std::size_t i=0; i<result.size(); ++i)
        to_value[i] = result[i];

      cell_setter(cell, to_value);

    } // dual_box_flux_to_cell


    /** @brief Interpolates normal components of the flux defined on each facet to all cells in the mesh (or segment).
     *
     * @param device           The ViennaSHE device on which to interpolate.
     * @param cell_setter      Functor for storing the interpolated flux vector
     * @param facet_accessor   Functor for accessing the normal components of the flux
     */
    template <typename DeviceType, typename CellSetter, typename FacetAccessor>
    void dual_box_flux_to_cell(DeviceType const & device, CellSetter & cell_setter, FacetAccessor const & facet_accessor)
    {
      /*
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennagrid::result_of::cell<MeshType>::type                   CellType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type       CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type          CellIterator;

      typedef typename viennagrid::result_of::const_facet_range<CellType>::type      FacetOnCellContainer;

      CellContainer cells(device.mesh());
      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        FacetOnCellContainer facets(*cit);
        dual_box_flux_to_cell(*cit, facets, cell_setter, facet_accessor);
      }*/

      throw std::runtime_error("dual_box_flux_to_cell(): Not implemented!");
    }

  } // util
} // viennashe


#endif /* VIENNASHE_UTIL_TRANSLATE_QUANTITY_HPP */

