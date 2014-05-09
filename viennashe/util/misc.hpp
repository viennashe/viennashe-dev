#ifndef VIENNASHE_UTIL_MISC_HPP
#define VIENNASHE_UTIL_MISC_HPP

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
#include <iomanip>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cmath>

// viennagrid
#include "viennagrid/forwards.hpp"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/algorithm/interface.hpp"
#include "viennagrid/algorithm/boundary.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/util/exception.hpp"
#include "viennashe/exception.hpp"

#include "viennashe/materials/all.hpp"

/** @file viennashe/util/misc.hpp
    @brief Miscellaneous utilities
*/

namespace viennashe
{
  /** @brief A collection of utilities used within ViennaSHE. */
  namespace util
  {

    //
    // Helper
    //

    /**
     * @brief A functor which holds a single value. Used in most of the postprocessing routines.
     * @tparam ValueType The type of the value to hold, typically a double.
     */
    template <typename ValueType>
    class value_holder_functor
    {
    public:
      typedef ValueType   value_type;

      template <typename T>
      void operator()(T const &, value_type val) const { value_ = val; }

      value_type operator()() const { return value_; }

    private:
      mutable value_type value_;
    };

    //
    // Averaging
    //

    /** @brief A functor which computes the arithmetic average of all entries in a container.
     *
     * The container is required to provide STL-compatible forward-iterators and a member function size() returning the number of elements.
     */
    struct arithmetic_averaging
    {
      template <typename ContainerType>
      double operator()(ContainerType const & cont) const
      {
        double ret = 0;
        for (typename ContainerType::const_iterator cit  = cont.begin();
                                                    cit != cont.end();
                                                  ++cit)
              ret += *cit / cont.size();

        return ret;
      }
    };

    /** @brief A functor which computes the geometric average of all entries in a container.
     *
     * The container is required to provide STL-compatible forward-iterators and a member function size() returning the number of elements.
     * Note that all entries need to be non-negative for even container sizes.
     */
    struct geometric_averaging
    {
      template <typename ContainerType>
      double operator()(ContainerType const & cont) const
      {
        double ret = 1.0;
        for (typename ContainerType::const_iterator cit  = cont.begin();
                                                    cit != cont.end();
                                                  ++cit)
        {
          ret *= std::pow(*cit, 1.0 / cont.size());
        }

        return ret;
      }
    };

    /** @brief A functor which returns the result of a logic-or chain consisting of all elements of the supplied container.
     *
     * The container is required to provide STL-compatible forward-iterators and a member function size() returning the number of elements.
     */
    struct logic_or_averaging
    {
      template <typename ContainerType>
      bool operator()(ContainerType const & cont) const
      {
        bool ret = false;
        for (typename ContainerType::const_iterator cit  = cont.begin();
                                                    cit != cont.end();
                                                  ++cit)
        {
          ret |= *cit;
        }

        return ret;
      }
    };


    //
    // Filtered setter
    //
    /** @brief A compound functor which applies the provided functor only to the argument if the supplied filter evaluates to true
     *
     *  @tparam FilterType      Type of the filter functor. Needs to return 'true' if the action functor should be applied to the element
     *  @tparam FunctorType     The action functor operating on the supplied argument.
     */
    template <typename FilterType, typename FunctorType>
    class filtered_functor
    {
      public:
        filtered_functor(FilterType const & fil, FunctorType & fun) : filter_(fil), fun_(fun) {}

        template <typename T>
        void operator()(T const & t)
        {
          if (filter_(t))
            fun_(t);
        }

        template <typename T>
        void operator()(T & t)
        {
          if (filter_(t))
            fun_(t);
        }

      private:
        FilterType const & filter_;
        FunctorType & fun_;
    };

    /** @brief Convenience routine for creating a filtered functor out of the provided filter- and action-functor. Also see class filtered_functor. */
    template <typename FilterType, typename FunctorType>
    filtered_functor<FilterType, FunctorType>
    make_filtered_functor(FilterType const & filter, FunctorType & fun)
    {
      return filtered_functor<FilterType, FunctorType>(filter, fun);
    }


    /** @brief Namespace containing implementation details for functionality from viennashe::util. Typically not of interest for a library user. */
    namespace detail
    {
      /** @brief A functor assigning unity weight to the element provided. */
      struct unity_weight
      {
        template <typename U>
        double operator()(U const &) const { return 1.0; }
      };

    } //namespace detail


    /** @brief Helper function returning a const-pointer to the 'second cell' of a facet, or NULL if there is no second facet. The 'first' cell is passed to the function. */
    template <typename MeshT, typename FacetT, typename CellT>
    CellT const * get_other_cell_of_facet(MeshT const & mesh, FacetT const & facet, CellT const & cell)
    {
      typedef typename viennagrid::result_of::const_coboundary_range<MeshT, FacetT, CellT>::type    CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                  CellOnFacetIterator;

      CellOnFacetContainer cells_on_facet(mesh, viennagrid::handle(mesh, facet));
      CellOnFacetIterator  cofit = cells_on_facet.begin();

      if ( &(*cofit) == &(cell))  //one of the two vertices of the edge is different from 'cell'
        ++cofit;

      if (cofit == cells_on_facet.end()) // Only one cell attached to facet
        return NULL;

      return &(*cofit);
    }


    /** @brief Helper function returning a const-pointer to the 'second cell' of a facet, or NULL if there is no second facet. The 'first' cell is passed to the function. */
    template <typename DeviceT, typename CellT>
    CellT const * get_connected_semiconductor_cell(DeviceT const & device, CellT const & cell)
    {
      typedef typename viennagrid::result_of::const_facet_range<CellT>::type        FacetOnCellContainer;
      typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

      FacetOnCellContainer facets_on_cell(cell);
      for (FacetOnCellIterator focit = facets_on_cell.begin();
          focit != facets_on_cell.end();
          ++focit)
      {
        // check whether cell is adjacent to a semiconductor cell:
        CellT const * other_cell = viennashe::util::get_other_cell_of_facet(device.mesh(), *focit, cell);
        if (!other_cell)
          continue;

        if (viennashe::materials::is_semiconductor(device.get_material(*other_cell)))
          return other_cell;
      }

      return NULL;
    }


    ///////////////////// Others ////////////////////////

    // empty default template
    template <bool b>
    struct cpp03_static_assert {};

    // template specialized on true
    template <>
    struct cpp03_static_assert<true>
    {
      static void check() {}
    };

    template <typename T, typename U>
    struct is_equal
    {
      static const bool value = false;
    };

    template <typename T>
    struct is_equal<T, T>
    {
      static const bool value = true;
    };


    /** @brief A functor-style wrapper for a spatial quantity (typically potential, electron density or hole density). Typically evaluated on each vertex only.
     *
     * @tparam DeviceType       Type of the underlying device
     * @tparam dim              Topological dimension of the elements on which the quantity is defined (0: vertices)
     * @tparam ValueType        Value type of the quantity (usually double or bool)
     */
    template <typename DeviceType, typename ElementTagT, typename ValueType = double>
    class spatial_quantity_wrapper
    {
        typedef spatial_quantity_wrapper self_type;

        typedef typename DeviceType::mesh_type       MeshType;

      public:
        typedef typename viennagrid::result_of::element<MeshType, ElementTagT>::type       element_type;
        typedef ValueType      value_type;

        /** @brief Copy over all values from the accessors so that this object can be easily passed around */
        template <typename VectorType, typename IndexKeyType, typename BoundaryValueAccessor>
        spatial_quantity_wrapper(DeviceType const & device,
                                 VectorType const & quantity_values,
                                 IndexKeyType const & index_array,
                                 BoundaryValueAccessor const & bnd_accessor)
          : values_(viennagrid::elements<ElementTagT>(device.mesh()).size()), vec_(quantity_values.size())
        {
          cpp03_static_assert< is_equal<value_type,
                                        typename VectorType::value_type>::value
                             >::check();

          typedef typename viennagrid::result_of::const_element_range<MeshType, ElementTagT>::type     ElementContainer;
          typedef typename viennagrid::result_of::iterator<ElementContainer>::type                       ElementIterator;
          typedef typename element_type::id_type     id_type;

          // Iterate over all elements and copy values over
          ElementContainer elements(device.mesh());
          for (ElementIterator it = elements.begin(); it != elements.end(); ++it)
          {
            id_type element_id = it->id();
            long index = index_array.at(element_id);

            if (index >= 0)
              values_.at(element_id) = quantity_values[index];
            else
              values_.at(element_id) = bnd_accessor(*it);
          }

          for (std::size_t i=0; i<quantity_values.size(); ++i)
            vec_[i] = quantity_values[i];
        }

        /** @brief The functor interface. */
        value_type operator()(element_type const & t) const
        {
          return values_.at(t.id());
        }

        /** @brief Returns the internal vector with the quantities. Allows to externally adjust the quantity. */
        std::vector<double> const & vector() const { return vec_; }

      private:
        std::vector<value_type> values_;
        std::vector<double> vec_;
    };

    /** @brief A functor-style wrapper for a spatial quantity which is externally prescribed by the user.
     *
     * @tparam DeviceType       Type of the underlying device
     * @tparam dim              Topological dimension of the elements on which the quantity is defined (0: vertices)
     * @tparam ValueType        Value type of the quantity (usually double or bool)
     */
    template <typename DeviceType, typename ElementTagT, typename ValueType = double>
    class spatial_quantity
    {
        typedef spatial_quantity self_type;

        typedef typename DeviceType::mesh_type       MeshType;

      public:
        typedef typename viennagrid::result_of::element<MeshType, ElementTagT>::type       element_type;
        typedef ValueType      value_type;

        /** @brief Copy over all values from the accessors so that this object can be easily passed around */
        spatial_quantity(DeviceType const & device, ValueType default_value = ValueType()) : values_(viennagrid::cells(device.mesh()).size(), default_value) {}

        /** @brief The functor interface. */
        value_type operator()(element_type const & t) const
        {
          return values_.at(t.id());
        }

        void set(element_type const & t, ValueType val)
        {
          values_.at(t.id()) = val;
        }

      private:
        std::vector<value_type> values_;
    };


    /**
     * @brief Returns the formatted number as string with the given length, where spaces are used as padding character.
     * @param number The number to be formatted
     * @param width Non negative number
     * @return Formatted number as string
     */
    template <typename T>
    std::string format_number_to_width(T number, int width)
    {
        std::stringstream s;
        s << std::setfill(' ') << std::setw(width) << number;
        return s.str();
    }

    /** @brief Computes the infimum norm of a subvector
     *
     * @param v             The vector
     * @param index_start   The start index (inclusive)
     * @param index_stop    The stop index (exclusive)
     */
    template <typename VectorType>
    typename VectorType::value_type norm_inf(VectorType const & v, std::size_t index_start, std::size_t index_stop)
    {
      typename VectorType::value_type ret = 0;
      for (std::size_t i=index_start; i<index_stop; ++i)
        ret = std::max(std::abs(v[i]), ret);

      return ret;
    }


  } //namespace util
} //namespace viennashe
#endif
