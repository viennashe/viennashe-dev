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
#include "viennagrid/viennagrid.h"

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
      viennagrid_mesh mesh = device.mesh();

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
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
          viennagrid_element_id *facets_on_cell_begin, *facets_on_cell_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(mesh, *cit, cell_dim - 1, &facets_on_cell_begin, &facets_on_cell_end));
          for (viennagrid_element_id *focit  = facets_on_cell_begin;
                                      focit != facets_on_cell_end;
                                    ++focit)
          {
            viennagrid_element_id *other_cell;
            viennashe::util::get_other_cell_of_facet(mesh, *focit, *cit, &other_cell);
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

            for (viennagrid_element_id *focit  = facets_on_cell_begin;
                                        focit != facets_on_cell_end;
                                      ++focit)
            {
              viennagrid_element_id *other_cell;
              viennashe::util::get_other_cell_of_facet(mesh, *focit, *cit, &other_cell);
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
      viennagrid_mesh mesh = device.mesh();

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

      //
      // Ensure that expansion order between vertices does not vary too much. If difference is too large, reduce order of 'higher' vertex.
      //
      for (long i=0; i<10; ++i)
      {
        //log::info<log_smooth_expansion_order>() << "* smooth_expansion_order(): Vertex expansion order smoothing, run " << i+1 << "..." << std::endl;
        bool something_changed = false;

        viennagrid_element_id *facets_begin, *facets_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim-1, &facets_begin, &facets_end));
        for (viennagrid_element_id *fit  = facets_begin;
                                    fit != facets_end;
                                  ++fit)
        {
          viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

          viennagrid_element_id c1 = cells_on_facet_begin[0];

          if (cells_on_facet_begin + 1 == cells_on_facet_end)
            continue;

          viennagrid_element_id c2 = cells_on_facet_begin[1];

          for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if ( quan.get_expansion_order(c1, index_H) > quan.get_expansion_order(c2, index_H) + 2)
            {
              something_changed = true;
              quan.set_expansion_order(c1, index_H, quan.get_expansion_order(c2, index_H) + 2);
            }
            if ( quan.get_expansion_order(c2, index_H) > quan.get_expansion_order(c1, index_H) + 2)
            {
              something_changed = true;
              quan.set_expansion_order(c2, index_H, quan.get_expansion_order(c1, index_H) + 2);
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
      viennagrid_mesh mesh = device.mesh();

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

      viennagrid_element_id *facets_begin, *facets_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim - 1, &facets_begin, &facets_end));
      for (viennagrid_element_id *fit = facets_begin;
                                  fit != facets_end;
                                ++fit)
      {
        viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

        for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
        {
          std::size_t order = quan.get_expansion_order(cells_on_facet_begin[0], index_H);

          if (cells_on_facet_begin + 1 != cells_on_facet_end)
            order = std::max<std::size_t>(order, quan.get_expansion_order(cells_on_facet_begin[1], index_H));
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
      template <typename DeviceType>
      void distribute_expansion_orders_impl(DeviceType const & device,
                                            viennashe::she::unknown_she_quantity<double> & quan,
                                            viennashe::config const & conf)
      {
        viennagrid_mesh mesh = device.mesh();

        std::size_t Lmax = static_cast<std::size_t>(conf.max_expansion_order());

        //
        // write expansion order on cells (even unknowns)
        //

        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

        viennagrid_element_id *cells_begin, *cells_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
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
        viennagrid_element_id *facets_begin, *facets_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim - 1, &facets_begin, &facets_end));
        for (viennagrid_element_id *fit = facets_begin;
                                    fit != facets_end;
                                  ++fit)
        {
          viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if ( conf.adaptive_expansions() )
            {
              std::size_t order = quan.get_expansion_order(cells_on_facet_begin[0], index_H);

              if (cells_on_facet_begin + 1 != cells_on_facet_end)
                order = std::max<std::size_t>(order, quan.get_expansion_order(cells_on_facet_begin[1], index_H));
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
    template <typename DeviceType>
    void distribute_expansion_orders(DeviceType & device,
                                     viennashe::she::unknown_she_quantity<double> & quan,
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
