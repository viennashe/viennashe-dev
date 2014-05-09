#ifndef VIENNASHE_SHE_EXPANSION_ORDER_HPP
#define VIENNASHE_SHE_EXPANSION_ORDER_HPP

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

// std
#include <iostream>

// viennagrid
#include "viennagrid/element/element.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/mesh/coboundary_iteration.hpp"

// viennashe
#include "viennashe/util/misc.hpp"
#include "viennashe/forwards.h"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/config.hpp"
#include "viennashe/she/assemble_common.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/she/timestep_quantities.hpp"
#include "viennashe/she/she_quantity.hpp"

#include "viennashe/she/exception.hpp"

#include "viennashe/util/filter.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"

/** @file viennashe/she/expansion_order.hpp
    @brief Provides all routines for distributing unknown indices (dofs) over the device. Includes expansion order smoothing routines.
*/

namespace viennashe
{
  namespace she
  {


    /** @brief Assigns the various expansion orders to the device. Uniform expansion order implemented here
     *
     * @param device        The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan          The SHE quantities
     */
    template <typename DeviceType, typename SHEQuantity>
    void lower_expansion_order_near_bandedge(DeviceType const & device, SHEQuantity & quan)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

      typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
      typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

      MeshType const & mesh = device.mesh();

      CellContainer cells(mesh);
      for (CellIterator cit = cells.begin();
          cit != cells.end();
          ++cit)
      {

        for (std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          if (quan.get_kinetic_energy(*cit, index_H) <= 0)   // In forbidden band there is expansion order 1
          {
            quan.set_expansion_order(*cit, index_H, 1);
            continue;
          }


          //
          // Make sure that there is no vertex close (w.r.t. energy and space) to the band edge with order above 1
          //

          bool all_neighbors_have_positive_kinetic_energy = true;
          //bring all neighbors of the vertex to order 1
          FacetOnCellContainer facets_on_cell(*cit);
          for (FacetOnCellIterator focit = facets_on_cell.begin();
              focit != facets_on_cell.end();
              ++focit)
          {
            CellType const *other_cell = viennashe::util::get_other_cell_of_facet(mesh, *focit, *cit);
            if (!other_cell)
              continue;

            if (quan.get_kinetic_energy(*other_cell, index_H-1) <= 0)   // In forbidden band there is expansion order 1
              all_neighbors_have_positive_kinetic_energy = false;
          }


          if (all_neighbors_have_positive_kinetic_energy)
            break;
          else
          {
            //set expansion order of all neighbors to one:
            quan.set_expansion_order(*cit, index_H, 1);

            for (FacetOnCellIterator focit = facets_on_cell.begin();
                focit != facets_on_cell.end();
                ++focit)
            {
              CellType const *other_cell = viennashe::util::get_other_cell_of_facet(mesh, *focit, *cit);
              if (!other_cell)
                continue;

              quan.set_expansion_order(*other_cell, index_H, 1);
            }
          }
        }
      }
    }



    /** @brief Smoothes the expansion order over the device such that neighboring vertices differ in their expansion order by less than two.
     *
     * @param device        The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan          The SHE quantities
     */
    template <typename DeviceType,
              typename SHEQuantity>
    void smooth_expansion_order(DeviceType const & device,
                                SHEQuantity  & quan)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;
      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;


      MeshType const & mesh = device.mesh();

      //
      // Ensure that expansion order between vertices does not vary too much. If difference is too large, reduce order of 'higher' vertex.
      //
      for (long i=0; i<10; ++i)
      {
        //log::info<log_smooth_expansion_order>() << "* smooth_expansion_order(): Vertex expansion order smoothing, run " << i+1 << "..." << std::endl;
        bool something_changed = false;

        FacetContainer facets(mesh);
        for (FacetIterator fit = facets.begin();
             fit != facets.end();
             ++fit)
        {
          CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), fit.handle()));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          CellType const & c1 = *cofit;
          ++cofit;

          if (cofit == cells_on_facet.end())
            continue;

          CellType const *c2_ptr = &(*cofit);

          for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if ( quan.get_expansion_order(c1, index_H) > quan.get_expansion_order(*c2_ptr, index_H) + 2)
            {
              something_changed = true;
              quan.set_expansion_order(c1, index_H, quan.get_expansion_order(*c2_ptr, index_H) + 2);
            }
            if ( quan.get_expansion_order(*c2_ptr, index_H) > quan.get_expansion_order(*c2_ptr, index_H) + 2)
            {
              something_changed = true;
              quan.set_expansion_order(*c2_ptr, index_H, quan.get_expansion_order(*c2_ptr, index_H) + 2);
            }
          } //for index_H
        } //for eit

        if (!something_changed)
          break;
      }
    }

    /** @brief Writes expansion orders on edges. Edge expansion order is given as max(order(v1), order(v2)), where v1 and v2 denote the cells of the facet
     *
     * @param device        The device on which simulation is carried out
     * @param quan          The SHE quantities
     */
    template <typename DeviceType,
              typename SHEQuantity>
    void write_expansion_order_to_facets(DeviceType & device,
                                        SHEQuantity  & quan)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;
      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

      MeshType const & mesh = device.mesh();

      FacetContainer facets(mesh);
      for (FacetIterator fit = facets.begin();
           fit != facets.end();
           ++fit)
      {
        CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), fit.handle()));

        CellOnFacetIterator cofit = cells_on_facet.begin();
        CellType const & c1 = *cofit;
        ++cofit;
        CellType const *c2_ptr = NULL;

        if (cofit != cells_on_facet.end())
          c2_ptr = &(*cofit);

        for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
        {
          std::size_t order = quan.get_expansion_order(c1, index_H);

          if (c2_ptr != NULL)
            order = std::max<std::size_t>(order, quan.get_expansion_order(*c2_ptr, index_H));
          quan.set_expansion_order(*fit, index_H, order);
        }
      }

    }


    namespace detail
    {

      /** @brief Assigns the various expansion orders to the device. This is the implementation.
       *
       * @param device        The device (includes a ViennaGrid mesh) on which simulation is carried out
       * @param quan          The unkown SHE quantities on edges and vertices
       * @param conf          The simulator configuration
      */
      template <typename DeviceType, typename VertexT, typename EdgeT>
      void distribute_expansion_orders_impl(DeviceType const & device,
                                            viennashe::she::unknown_she_quantity<VertexT, EdgeT> & quan,
                                            viennashe::config const & conf)
      {
        typedef typename DeviceType::mesh_type              MeshType;

        typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;
        typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;

        typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
        typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

        typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

        typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
        typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

        MeshType const & mesh = device.mesh();

        std::size_t Lmax = static_cast<std::size_t>(conf.max_expansion_order());

        //
        // write expansion order on cells (even unknowns)
        //
        CellContainer cells(mesh);
        for (CellIterator cit = cells.begin();
            cit != cells.end();
            ++cit)
        {
          if ( (conf.adaptive_expansions() == false)
              || (Lmax == 1) )
          {
            for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
              quan.set_expansion_order(*cit, index_H, Lmax);
          }
        }

        //
        // Write expansion order on facets (odd unknowns)
        //
        FacetContainer facets(mesh);
        for (FacetIterator fit = facets.begin();
             fit != facets.end();
             ++fit)
        {
          CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), fit.handle()));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          CellType const & c1 = *cofit;
          ++cofit;
          CellType const *c2_ptr = NULL;

          if (cofit != cells_on_facet.end())
            c2_ptr = &(*cofit);

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if ( conf.adaptive_expansions() )
            {
              std::size_t order = quan.get_expansion_order(c1, index_H);

              if (c2_ptr != NULL)
                order = std::max<std::size_t>(order, quan.get_expansion_order(*c2_ptr, index_H));
              quan.set_expansion_order(*fit, index_H, order);
            }
            else
              quan.set_expansion_order(*fit, index_H, Lmax);
          }
        }

      } //distribute_expansion_orders_impl




    } //namespace detail


    /** @brief Assigns the various expansion orders to the device. This is the public interface.
     *
     * @param device        The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan          The SHE unkown quantities on edges and vertices
     * @param conf          The simulator configuration
     */
    template <typename DeviceType,
              typename VertexT, typename EdgeT>
    void distribute_expansion_orders(DeviceType & device,
                                     viennashe::she::unknown_she_quantity<VertexT, EdgeT> & quan,
                                     viennashe::config const & conf)
    {
      detail::distribute_expansion_orders_impl(device, quan, conf);

      //
      // At band edge use first order expansions
      //
      lower_expansion_order_near_bandedge(device, quan);

      //
      // Bring edges in consistent state:
      //
      smooth_expansion_order(device, quan);
    }

    template <typename DeviceType>
    void distribute_expansion_orders(DeviceType const & device,
                                     timestep_quantities<DeviceType> & quantities,
                                     viennashe::config const & conf)
    {
      typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_she_quantity_type  SHEUnknownType;

      SHEUnknownType & f_n = quantities.electron_distribution_function();
      SHEUnknownType & f_p = quantities.hole_distribution_function();

      // Setup energies for electrons:
      if (conf.with_electrons() && conf.get_electron_equation() == EQUATION_SHE)
        distribute_expansion_orders(device, f_n, conf);

      // Setup energies for holes:
      if (conf.with_holes() && conf.get_hole_equation() == EQUATION_SHE)
        distribute_expansion_orders(device, f_p, conf);
    }


  }
}

#endif
