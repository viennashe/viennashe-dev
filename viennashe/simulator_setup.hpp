#ifndef VIENNASHE_SIMULATOR_SETUP_HPP
#define VIENNASHE_SIMULATOR_SETUP_HPP

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

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <limits>

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/accessors.hpp"
#include "viennashe/setters.hpp"
#include "viennashe/device.hpp"
#include "viennashe/exception.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/log/log.hpp"

#include "viennagrid/viennagrid.h"


/** @file  viennashe/simulator_setup.hpp
    @brief Defines a set of setup routines and checkers for simulators.
*/

namespace viennashe
{
  namespace detail
  {

    /** @brief Ensures that a quantity is available on n-cells of a certain dimension. If this is not the case, the quantity is transferred (averaged/interpolated) from cells of other dimension. An exception is thrown if the fallback-transfer fails.
     *
     * @tparam dim            Topological dimension of the n-cell for which the quantity should be available
     * @tparam DeviceType     Type of the device
     * @tparam AccessorType   An accessor for the quantity for an n-cell
     * @tparam CheckerType    Type of the checker functor, which checks for a valid range
     * @tparam FilterType     A filter for the n-cells to actually consider
     */
    template<typename DeviceType, typename AccessorType, typename FilterType, typename CheckerType>
    bool check_quantity_on_ncell(DeviceType const & device,
                                 viennagrid_dimension topo_dim,
                                 AccessorType const & accessor,
                                 FilterType const & filter,
                                 CheckerType const & checker)
    {
      typedef typename DeviceType::mesh_type           MeshType;

      viennagrid_element_id *elements_begin, *elements_end;
      viennagrid_mesh_elements_get(device.mesh(), topo_dim, &elements_begin, &elements_end);

      bool quantity_okay = true;

      // Check whether quantities are okay:
      for (viennagrid_element_id *it  = elements_begin;
                                  it != elements_end;
                                ++it)
      {
        if (!filter(*it))
          continue;

        if (!checker(accessor(*it)))
        {
          quantity_okay = false;
          break;
        }
      }

      return quantity_okay;
    }


    /** @brief Ensures that a quantity is available on n-cells of a certain dimension. If this is not the case, the quantity is transferred (averaged/interpolated) from cells of other dimension.
     *
     * @tparam TagT                  Topological tag of the n-cell for which the quantity should be available
     * @tparam FallbackTagT          Topological tag of the n-cells from which the data should be taken otherwise
     * @tparam DeviceType            Type of the device
     * @tparam AccessorType          An accessor for the quantity for an n-cell
     * @tparam CheckerType           Type of the checker functor, which checks for a valid range
     * @tparam FallbackAccessorType  Type of the accessor to be used if the quantity is not set up properly on the mesh
     * @tparam FallbackSetterType    Type of the setter to be used if the quantity is not set up properly on the mesh
     * @tparam FallbackAveragerType  Type of the averager to be used if the quantity is not set up properly on the mesh
     * @tparam FallbackFilterType    Type of the filter to be used for grabbing data from fallback cells. Allows to exclude specific fallback cells
     * @tparam FallbackCheckerType   Type of the checker to ensure that the quantity is set up correctly after the fallback
     *
     * @return True if a fallback has been used, otherwise false.
     */
    template <typename TagT, typename FallbackTagT, typename DeviceType,
              typename AccessorType, typename FilterType, typename CheckerType,
              typename FallbackAccessorType, typename FallbackSetterType, typename FallbackAveragerType, typename FallbackFilterType, typename FallbackCheckerType>
    bool setup_and_check_quantity_on_ncell_with_fallback(DeviceType & device,
                                                         AccessorType const & accessor,
                                                         FilterType const & filter,
                                                         CheckerType const & checker,
                                                         FallbackAccessorType const & fallback_accessor,
                                                         FallbackSetterType         & fallback_setter,
                                                         FallbackAveragerType const & fallback_averager,
                                                         FallbackFilterType   const & fallback_filter,
                                                         FallbackCheckerType  const & fallback_checker
                                                        )
    {
      bool quantity_okay = check_quantity_on_ncell<TagT>(device, accessor, filter, checker);

      // if quantity is not yet okay, use fallback
      if (!quantity_okay)
      {
        // transfer quantity
        /*viennagrid::quantity_transfer<FallbackTagT, TagT>(device.mesh(),
                                                          fallback_accessor, fallback_setter, // Accessor (in) and Setter (out)
                                                          fallback_averager,           // Averaging
                                                          fallback_filter, filter);    // Filters (in, out)

        //check again:
        check_quantity_on_ncell<TagT>(device, accessor, filter, fallback_checker);
        return true;*/

        throw std::runtime_error("setup_and_check_quantity_on_ncell_with_fallback(): Fallback not implemented!");
      }

      return false;
    }

    /** @brief Helper routine for creating the exception thrown if the donator doping is found to be invalid on the device */
    inline invalid_doping_in_device_exception get_invalid_doping_exception(viennashe::carrier_type_id ctype)
    {
      if (ctype == viennashe::ELECTRON_TYPE_ID)
        return invalid_doping_in_device_exception("Check failed: Invalid donator doping in device!");
      else if (ctype == viennashe::HOLE_TYPE_ID)
        return invalid_doping_in_device_exception("Check failed: Invalid acceptor doping in device!");
      else
        return invalid_doping_in_device_exception("Check failed: Invalid UNKOWN doping in device!");
    }

  } // namespace detail


  //
  // Doping
  //

  /** @brief Sets up donator or acceptor doping on vertices. If doping is provided on cells, it is interpolated accordingly.
   * @tparam DeviceType The device type (cf. device.hpp)
   * @param device A non-const reference to the device
   * @param ctype The carrier type identifier (cf. forwards.h)
   */
  template <typename DeviceType>
  void setup_doping_on_vertices(DeviceType & device, viennashe::carrier_type_id ctype)
  {
    /*
    typedef typename viennagrid::result_of::cell_tag<typename DeviceType::mesh_type>::type  CellTag;

    // fallback functors:
    doping_setter<DeviceType> fallback_setter(device, ctype);
    util::geometric_averaging fallback_averager;
    util::semiconductor_filter<DeviceType> fallback_filter(device);

    // standard check functors:
    doping_accessor<DeviceType> doping_accessor(device, ctype);
    util::vicinity_filter<DeviceType,
                          CellTag,
                          util::semiconductor_filter<DeviceType> > filter(device, fallback_filter);
    util::range_checker<double> checker(1.0, 1e32);

    bool fallback_used = detail::setup_and_check_quantity_on_ncell_with_fallback<viennagrid::vertex_tag, CellTag>(
                             device, doping_accessor, filter, checker,
                             doping_accessor, fallback_setter, fallback_averager, fallback_filter,
                             util::make_checker_with_exception(checker, detail::get_invalid_doping_exception(ctype))
                            );
    if (fallback_used)
    {
      log::info() << "* setup_doping_on_vertices(): Doping defined on cells successfully transferred to vertices. " << std::endl;
    }
    else
    {
      log::info() << "* setup_doping_on_vertices(): Doping defined on vertices passes consistency check." << std::endl;
    } */

    throw std::runtime_error("setup_doping_on_vertices(): Not implemented!");

  }

  /** @brief Sets up donator and acceptor doping on vertices. If doping is provided on cells, it is interpolated accordingly.
   * @tparam DeviceType The device type (cf. device.hpp)
   * @param device A non-const reference to the device
   */
  template <typename DeviceType>
  void setup_doping_on_vertices(DeviceType & device)
  {
    setup_doping_on_vertices(device, viennashe::ELECTRON_TYPE_ID);
    setup_doping_on_vertices(device, viennashe::HOLE_TYPE_ID);
  }


  /**
   * @brief Calculates the distances in each semiconductor segment to the next neighboring insulator segment on each vertex.
   *        Uses ViennaData to store the calculated distances
   *
   * @tparam DeviceType The device type (cf. device.hpp)
   * @param device A non-const reference to the device
   * @todo What about corners ? (think of an L-shaped insulator and a semiconductor in its corner)
   */
  /*
  template <typename DeviceType>
  void setup_insulator_distances(DeviceType & device)
  {
    typedef typename DeviceType::mesh_type                                          MeshType;
    typedef typename viennagrid::result_of::segmentation<MeshType>::type          SegmentationType;
    typedef typename viennagrid::result_of::segment_handle<SegmentationType>::type  SegmentType;
    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;

    typedef typename viennagrid::result_of::const_cell_range<SegmentType>::type       CellOnSegmentContainer;
    typedef typename viennagrid::result_of::const_vertex_range<SegmentType> ::type    VertexOnSegmentContainer;
    typedef typename viennagrid::result_of::iterator<VertexOnSegmentContainer>::type  VertexIterator;

    log::info() << "* setup_insulator_distances(): Calculating interface distances ..." << std::endl;

    // for every semiconductor segment ...
    for ( typename SegmentationType::const_iterator seg_it  = device.segmentation().begin();
                                                    seg_it != device.segmentation().end();
                                                  ++seg_it)
    {
      const SegmentType & segment = *seg_it;
      CellOnSegmentContainer cells(segment);

      if(cells.size() > 0 && viennashe::materials::is_semiconductor(device.get_material(cells[0]))) // exclude non semiconductors
      {
        // ... loop over every vertex in that segment ...
        VertexOnSegmentContainer vertices(segment);
        for (VertexIterator vit = vertices.begin(); vit != vertices.end(); vit++)
        {
          double mindist = std::numeric_limits<double>::max();
          // ... loop over all insulator segments ...
          for ( typename SegmentationType::const_iterator seg_it2  = device.segmentation().begin();
                                                          seg_it2 != device.segmentation().end();
                                                        ++seg_it2)
          {
            const SegmentType & insulator_segment = *seg_it2;
            CellOnSegmentContainer insulator_cells(insulator_segment);

            if(insulator_cells.size() == 0 || &segment == &insulator_segment) continue; // skip empty segments
            if(!viennashe::materials::is_insulator(device.get_material(insulator_cells[0]))) continue; // skip non-insulator segments

            // ... find the shortest distance from the vertex to the insulator segment and cache it ...
            throw "Fix me!";
            const PointType pt; // FIX = viennagrid::closest_points(*vit, insulator_segment).second - viennagrid::point(*vit);
            const double dist  = viennagrid::norm_2(viennagrid::point(*vit) - pt);
            // ... if it is smaller than the previous distance
            if(dist < mindist)
            {
              //viennashe::log::debug() << "DISTANCE vector from " << sid_semi << " to "<< sid << " for vertex " << vit->id() << " => (" << pt << ") " << std::endl;
              device.set_vector_to_next_insulator(*vit, pt); // TODO: What about corners ? (think of an L-shaped insulator and a semiconductor in its corner)
              mindist = dist;
            }
          } // for insulator segments
        } // for vertices of semiconductor segment
      } // if current segment is semiconductor
    } // for segments

  } */

} //namespace viennashe

#endif
