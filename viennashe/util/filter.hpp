#ifndef VIENNASHE_UTIL_FILTER_HPP
#define VIENNASHE_UTIL_FILTER_HPP

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

#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/enable_if.hpp"

/** @file viennashe/util/filter.hpp
    @brief Defines several filter functors for the device. A filter functor returns true if the supplied argument is accepted, false otherwise.
*/

namespace viennashe
{
  /** @brief A collection of utilities used within ViennaSHE. */
  namespace util
  {

    /** @brief A trivial filter, all objects are accepted */
    struct any_filter
    {
      template <typename T>
      bool operator()(T const &) const { return true; }
    };

    /** @brief A trivial filter, no objects shall be accepted */
    struct none_filter
    {
      template <typename T>
      bool operator()(T const &) const { return false; }
    };

    /** @brief  A filter accepting semiconductor cells only */
    template <typename DeviceType>
    class semiconductor_filter
    {
      public:
        semiconductor_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id cell) const { return viennashe::materials::is_semiconductor(device_.get_material(cell)); }

      private:
        DeviceType const & device_;
    };

    /** @brief A simple filter to find cells and vertices which have potential boundary conditions */
    template <typename DeviceType>
    class contact_filter
    {
      public:
        contact_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id element) const { return device_.has_contact_potential(element); }

      private:
        DeviceType const & device_;
    };

    /** @brief The inverse filter to the contact_filter */
    template <typename DeviceType>
    class no_contact_filter
    {
      public:
        no_contact_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id element) const { return !device_.has_contact_potential(element); }

      private:
        DeviceType const & device_;
    };


    /** @brief A filter used to filter cells, which are neither semiconductors (by material) nor contacts (by potential boundary condition) */
    template <typename DeviceType>
    class no_contact_and_semiconductor_filter
    {
      public:
        no_contact_and_semiconductor_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id element) const
        {
          return !device_.has_contact_potential(element) && viennashe::materials::is_semiconductor(device_.get_material(element));
        }

      private:
        DeviceType const & device_;
    };



    namespace detail
    {
      template <long a, long b>
      struct is_lesser
      {
        enum { value = (a < b) };
      };
    } // detail

    /** @brief  A filter returning true if any of the neighboring ncells of dimension 'dim' evaluate to true. */     //TODO: Think about adding this generic filter to ViennaGrid
    /* TODO: Migrate to ViennaGrid 3.0
    template <typename DeviceType, typename ElementTagT, typename CheckerType>
    class vicinity_filter
    {
        typedef typename DeviceType::mesh_type     MeshType;

      public:
        vicinity_filter(DeviceType const & d, CheckerType const & c) : device_(d), checker_(c) {}

        // checks on ncell coboundary
        template <typename NCellType>
        typename viennashe::enable_if< detail::is_lesser<NCellType::tag::dim, ElementTagT::dim>::value, //topological dimension
                                       bool>::type
        operator()(NCellType const & cell) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<MeshType, NCellType, ElementTagT>::type    VicinityContainer;
          typedef typename viennagrid::result_of::iterator<VicinityContainer>::type                             VicinityIterator;

          VicinityContainer vicinity(device_.mesh(), viennagrid::handle(device_.mesh(), cell));
          for (VicinityIterator vit  = vicinity.begin();
                                vit != vicinity.end();
                              ++vit)
          {
            if (checker_(*vit))
              return true;
          }

          return false;
        }

        // checks on ncell boundary
        template <typename NCellType>
        typename viennashe::enable_if< (NCellType::tag::dim > ElementTagT::dim), //topological dimension
                                        bool>::type
        operator()(NCellType const & cell) const
        {
          typedef typename viennagrid::result_of::const_element_range<NCellType, ElementTagT>::type  VicinityContainer;
          typedef typename viennagrid::result_of::iterator<VicinityContainer>::type                  VicinityIterator;

          VicinityContainer vicinity = viennagrid::elements<ElementTagT>(cell);
          for (VicinityIterator vit  = vicinity.begin();
                                vit != vicinity.end();
                              ++vit)
          {
            if (checker_(*vit))
              return true;
          }

          return false;
        }

      private:
        DeviceType const & device_;
        CheckerType const & checker_;
    };*/


    /** @brief  A filter accepting metal cells only */
    template <typename DeviceType>
    class metal_filter
    {
      public:
        metal_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id element) const { return viennashe::materials::is_conductor(device_.get_material(element)); }

      private:
        DeviceType const & device_;
    };

    /** @brief  A filter accepting oxide cells only */
    template <typename DeviceType>
    class oxide_filter
    {
      public:
        oxide_filter(DeviceType const & d) : device_(d) {}

        bool operator()(viennagrid_element_id element) const { return viennashe::materials::is_insulator(device_.get_material(element)); }

      private:
        DeviceType const & device_;
    };



    // feel free to add further filters here


  } //namespace util
} //namespace viennashe
#endif
