#ifndef VIENNASHE_UTIL_POSTPROC_ELECTRIC_FIELD_HPP
#define VIENNASHE_UTIL_POSTPROC_ELECTRIC_FIELD_HPP

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
#include <fstream>
#include <vector>

// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/util/dual_box_flux.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"

/** @file viennashe/postproc/electric_field.hpp
    @brief Computes the electric field from a potential
 */

namespace viennashe
{
  namespace detail
  {
    /** @brief An accessor to the electric field along a given edge. Electrostatic potential required */

    template <typename DeviceType,
              typename PotentialAccessorType>
    struct electric_field_on_facet
    {
      typedef typename DeviceType::mesh_type MeshType;

      electric_field_on_facet(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) {}

      double operator()(viennagrid_element_id facet) const
      {
        viennagrid_dimension cell_dim;
        viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim);

        viennagrid_element_id *cells_begin, *cells_end;
        viennagrid_element_coboundary_elements(device_.mesh(), facet, cell_dim, &cells_begin, &cells_end);

        if (cells_end == cells_begin + 1) // facet on surface
          return 0;

        const double potential_center = potential_.get_value(cells_begin[0]);
        const double potential_outer  = potential_.get_value(cells_begin[1]);

        std::vector<viennagrid_numeric> centroid_0(3);
        std::vector<viennagrid_numeric> centroid_1(3);

        viennagrid_element_centroid(device_.mesh(), cells_begin[0], &(centroid_0[0]));
        viennagrid_element_centroid(device_.mesh(), cells_begin[1], &(centroid_1[0]));

        viennagrid_numeric distance;
        viennagrid_distance_2(cell_dim, &(centroid_0[0]), &(centroid_1[0]), &distance);

        const double Emag     = -(potential_outer - potential_center) / distance;

        return Emag;
      }
      private:
        DeviceType const & device_;
        PotentialAccessorType const & potential_;

    }; // electric_field_on_edge

  } // namespace detail

  /** @brief An accessor to the electric field on vertices and edges. Potential requiered */
  template < typename DeviceType, typename PotentialAccessorType >
  struct electric_field_wrapper
  {
  public:
    typedef std::vector<double> value_type;

    electric_field_wrapper(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) { }

    std::vector<double> operator()(viennagrid_element_id cell_or_facet) const
    {
      viennagrid_dimension cell_dim;
      viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim);

      viennagrid_dimension element_dim = viennagrid_topological_dimension_from_element_id(cell_or_facet);

      if (cell_dim == element_dim) // cell provided
      {
        typedef detail::electric_field_on_facet<DeviceType, PotentialAccessorType>  FieldOnFacetEvaluator;

        viennashe::materials::checker no_conductor_filter(MATERIAL_NO_CONDUCTOR_ID);

        std::vector<double> E(3);

        if (!no_conductor_filter(device_.get_material(cell_or_facet))) return E;

        viennashe::util::value_holder_functor<std::vector<double> > result;

        FieldOnFacetEvaluator facet_evaluator(device_, potential_);

        viennashe::util::dual_box_flux_to_cell(device_, cell_or_facet, result, facet_evaluator);

        return result();
      }
      else if (cell_dim == element_dim + 1) // facet provided
      {
        viennashe::detail::electric_field_on_facet<DeviceType, PotentialAccessorType> facet_eval(device_, potential_);
        std::vector<double> ret(3);
        ret[0] = facet_eval(cell_or_facet);
        return ret;
      }

      throw std::runtime_error("electric_field_wrapper::operator(): Invalid element dimension");

      std::vector<double> ret(3);
      return ret; //facet_eval(facet);
    }

    private:
      DeviceType const & device_;
      PotentialAccessorType const & potential_;
  }; // electric_field_wrapper


  /** @brief Convenience function for writing the electric field to a container
   *
   * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param potential        A quantity accessor to the potential
   * @param container        The container to be filled with values
   * @tparam ContainerType   A container type, for example: std::vector or std::deque
   */
  template <typename DeviceType,
            typename PotentialAccessor>
  void write_electric_field_to_quantity_field(DeviceType const & device,
                                              PotentialAccessor const & potential,
                                              viennagrid_quantity_field field)
  {
    electric_field_wrapper<DeviceType, PotentialAccessor> Efield(device, potential);

    viennashe::write_macroscopic_quantity_to_quantity_field(device, Efield, field);
  }

} // viennashe

#endif

