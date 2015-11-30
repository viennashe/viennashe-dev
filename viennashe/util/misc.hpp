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
#include <vector>

// viennagrid
#include "viennagrid/viennagrid.h"

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


    /** @brief Helper function returning a const-pointer to the 'second cell' of a facet, or NULL if there is no second cell. The 'first' cell is passed to the function. */
    template <typename MeshT>
    void get_other_cell_of_facet(MeshT const & mesh, viennagrid_element_id facet, viennagrid_element_id cell, viennagrid_element_id **other_cell)
    {
      viennagrid_dimension facet_dim = viennagrid_topological_dimension_from_element_id(facet);

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, facet, facet_dim + 1, &cells_begin, &cells_end));

      *other_cell = cells_begin;

      if (cell == cells_begin[0])
      {
        if (cells_end == cells_begin + 1)
          *other_cell = NULL;
        else
          *other_cell = cells_begin + 1;
      }
    }


    /** @brief Comment was wrong, needs to be rewritten */
    template <typename DeviceT, typename CellT>
    viennagrid_element_id * get_connected_semiconductor_cell(DeviceT const & device, CellT const & cell)
    {
      throw std::runtime_error("get_connected_semiconductor_cell(): Not implemented!");

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
    template <typename DeviceType, typename ValueType = double>
    class spatial_quantity_wrapper
    {
        typedef spatial_quantity_wrapper self_type;

        typedef typename DeviceType::mesh_type       MeshType;

      public:
        typedef ValueType      value_type;

        /** @brief Copy over all values from the accessors so that this object can be easily passed around */
        template <typename VectorType, typename IndexKeyType, typename BoundaryValueAccessor>
        spatial_quantity_wrapper(DeviceType const & device,
                                 VectorType const & quantity_values,
                                 IndexKeyType const & index_array,
                                 BoundaryValueAccessor const & bnd_accessor)
        {
          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

          viennagrid_int cell_count;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(device.mesh(), cell_dim, &cell_count));
          values_.resize(cell_count);


          // Iterate over all elements and copy values over
          viennagrid_element_id *cells_begin, *cells_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
          for (viennagrid_element_id *cit = cells_begin; cit != cells_end; ++cit)
          {
            viennagrid_element_id element_index = viennagrid_index_from_element_id(*cit);
            long index = index_array.at(element_index);

            if (index >= 0)
              values_.at(element_index) = quantity_values[index];
            else
              values_.at(element_index) = bnd_accessor(index);
          }
        }

        /** @brief The functor interface. */
        value_type operator()(viennagrid_element_id const & t) const
        {
          return values_.at(viennagrid_index_from_element_id(t));
        }

      private:
        std::vector<value_type> values_;
    };

    /** @brief A functor-style wrapper for a spatial quantity which is externally prescribed by the user.
     *
     * @tparam DeviceType       Type of the underlying device
     * @tparam dim              Topological dimension of the elements on which the quantity is defined (0: vertices)
     * @tparam ValueType        Value type of the quantity (usually double or bool)
     */
    template <typename DeviceType, typename ValueType = double>
    class spatial_quantity
    {
        typedef spatial_quantity self_type;

        typedef typename DeviceType::mesh_type       MeshType;

      public:
        typedef ValueType      value_type;

        /** @brief Copy over all values from the accessors so that this object can be easily passed around */
        spatial_quantity(DeviceType const & device, ValueType default_value = ValueType())
        {
          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

          viennagrid_int cell_count;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(device.mesh(), cell_dim, &cell_count));
          values_.resize(cell_count);
        }

        /** @brief The functor interface. */
        value_type operator()(viennagrid_element_id const & t) const
        {
          return values_.at(viennagrid_index_from_element_id(t));
        }

        void set(viennagrid_element_id const & t, ValueType val)
        {
          values_.at(viennagrid_index_from_element_id(t)) = val;
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
