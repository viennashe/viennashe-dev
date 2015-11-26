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
    /*
    namespace detail
    {
      // 1d case:
      template <typename WrappedConfigT>
      typename viennagrid::result_of::point< viennagrid::element<viennagrid::line_tag,WrappedConfigT> >::type
      outer_cell_normal_at_facet(viennagrid::element<viennagrid::line_tag,   WrappedConfigT> const & cell,
                                 viennagrid::element<viennagrid::vertex_tag, WrappedConfigT> const & facet)
      {
        typedef typename viennagrid::result_of::point< viennagrid::element<viennagrid::line_tag,WrappedConfigT> >::type  PointType;

        PointType centroid_cell  = viennagrid::centroid(cell);
        PointType centroid_facet = viennagrid::centroid(facet);

        if (centroid_cell[0] < centroid_facet[0])
          return PointType(1.0);
        return PointType(-1.0);
      }

      // 2d, triangles and quadrilaterals:
      template <typename CellTagT, typename WrappedConfigT>
      typename viennagrid::result_of::point< viennagrid::element<CellTagT, WrappedConfigT> >::type
      outer_cell_normal_at_facet(viennagrid::element<CellTagT,             WrappedConfigT> const & cell,
                                 viennagrid::element<viennagrid::line_tag, WrappedConfigT> const & facet)
      {
        typedef viennagrid::element<viennagrid::line_tag, WrappedConfigT>            FacetType;
        typedef typename viennagrid::result_of::point<FacetType>::type               PointType;
        typedef typename viennagrid::result_of::const_vertex_range<FacetType>::type  VertexRange;

        VertexRange facet_vertices(facet);
        PointType facet_vec = viennagrid::point(facet_vertices[1]) - viennagrid::point(facet_vertices[0]);
        PointType facet_centroid_vec = viennagrid::centroid(cell) - viennagrid::point(facet_vertices[0]);

        PointType edge_normal(facet_vec[1], -facet_vec[0]); // create one normal to line

        if (viennagrid::inner_prod(facet_centroid_vec, edge_normal) > 0) // normal vector is pointing into cell, hence flip
          edge_normal *= -1;

        return edge_normal / viennagrid::norm(edge_normal);
      }

      // 3d, implementation for both tetrahedra and hexahedra:
      template <typename CellTagT, typename FacetTagT, typename WrappedConfigT>
      typename viennagrid::result_of::point< viennagrid::element<CellTagT, WrappedConfigT> >::type
      outer_cell_normal_at_facet_3d(viennagrid::element<CellTagT,             WrappedConfigT> const & cell,
                                    viennagrid::element<FacetTagT, WrappedConfigT> const & facet)
      {
        typedef viennagrid::element<FacetTagT, WrappedConfigT>                       FacetType;
        typedef typename viennagrid::result_of::point<FacetType>::type               PointType;
        typedef typename viennagrid::result_of::const_vertex_range<FacetType>::type  VertexRange;

        VertexRange facet_vertices(facet);
        PointType facet_vec1 = viennagrid::point(facet_vertices[1]) - viennagrid::point(facet_vertices[0]);
        PointType facet_vec2 = viennagrid::point(facet_vertices[2]) - viennagrid::point(facet_vertices[0]);
        PointType facet_centroid_vec = viennagrid::centroid(cell) - viennagrid::point(facet_vertices[0]);

        PointType edge_normal = viennagrid::cross_prod(facet_vec1, facet_vec2); // create one normal to line

        if (viennagrid::inner_prod(facet_centroid_vec, edge_normal) > 0) // normal vector is pointing into cell, hence flip
          edge_normal *= -1;

        return edge_normal / viennagrid::norm(edge_normal);
      }

      // 3d, overload for tetrahedra
      template <typename CellTagT, typename WrappedConfigT>
      typename viennagrid::result_of::point< viennagrid::element<CellTagT, WrappedConfigT> >::type
      outer_cell_normal_at_facet(viennagrid::element<CellTagT,                 WrappedConfigT> const & cell,
                                 viennagrid::element<viennagrid::triangle_tag, WrappedConfigT> const & facet)
      {
        return outer_cell_normal_at_facet_3d(cell, facet);
      }

      // 3d, overload for hexahedra
      template <typename CellTagT, typename WrappedConfigT>
      typename viennagrid::result_of::point< viennagrid::element<CellTagT, WrappedConfigT> >::type
      outer_cell_normal_at_facet(viennagrid::element<CellTagT,                      WrappedConfigT> const & cell,
                                 viennagrid::element<viennagrid::quadrilateral_tag, WrappedConfigT> const & facet)
      {
        return outer_cell_normal_at_facet_3d(cell, facet);
      }


    } */

    /** @brief Returns the unit outer normal of a facet with respect to the provided cell */
    std::vector<double>
    outer_cell_normal_at_facet(viennagrid_mesh mesh,
                               viennagrid_element_id cell,
                               viennagrid_element_id facet)
    {
      std::vector<double> ret(3);

      throw std::runtime_error("outer_cell_normal_at_facet(): TODO: implement!");

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
    template <typename DeviceT, typename CellT, typename FacetContainerT, typename CellSetterT, typename FacetAccessorT>
    void dual_box_flux_to_cell(DeviceT     const & device,
                               CellT       const & cell,        FacetContainerT const & facets,
                               CellSetterT       & cell_setter, FacetAccessorT  const & facet_access)
    {
      /*
      typedef typename viennagrid::result_of::iterator<FacetContainerT>::type   FacetOnCellIterator;
      typedef typename DeviceT::mesh_type                                       MeshType;
      typedef typename FacetOnCellIterator::value_type                          FacetType;
      typedef typename viennagrid::result_of::point<CellT>::type                PointType;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellT>::type    CellOnFacetContainer;

      viennashe::math::dense_matrix<double> M(PointType::dim, PointType::dim);
      std::vector<double>                   b(PointType::dim);

      std::vector<PointType>  normals(facets.size());
      std::vector<double>     flux_contributions(facets.size());

      // for all edges connected with the current vertex
      std::size_t facet_ctr = 0;
      for (FacetOnCellIterator focit  = facets.begin();
                               focit != facets.end();
                             ++focit, ++facet_ctr )
      {
        flux_contributions[facet_ctr] = facet_access(*focit);
        normals[facet_ctr]            = outer_cell_normal_at_facet(cell, *focit); // vector along the edge pointing away from Vertex 'to'

        // flip orientation of flux contribution if global orientation is different:
        CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), *focit));
        if (&cells_on_facet[0] != &cell)
          flux_contributions[facet_ctr] *= -1.0;
      }

      // assemble mass matrix and rhs: M_i * V_i = rhs_i
      for ( std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i )
      {
        for ( std::size_t j = 0; j < static_cast<std::size_t>(PointType::dim); ++j )
        {
          for ( std::size_t k = 0; k < normals.size(); ++k )
            M(i, j) += normals[k][i] * normals[k][j];
        }
        for ( std::size_t k = 0; k < normals.size(); ++k )
          b[i] += normals[k][i] * flux_contributions[k];
      }

      std::vector<double> to_value(PointType::dim); // default: zero vector

      std::vector<double> result = viennashe::solvers::solve(M, b);
      for (std::size_t i=0; i<result.size(); ++i)
        to_value[i] = result[i];

      cell_setter(cell, to_value);*/

      throw std::runtime_error("dual_box_flux_to_cell(): Not implemented!");

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

