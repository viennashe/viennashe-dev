#ifndef VIENNASHE_UTIL_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP
#define VIENNASHE_UTIL_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP

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
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/checks.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/util/misc.hpp"
#include "viennashe/util/dual_box_flux.hpp"

/** @brief Computes the electric flux density from a potential
 * @file viennashe/postproc/electric_flux_density.hpp
 */

namespace viennashe
{
  namespace detail
  {
    /**
     * @brief Simple accessor to get the electric flux density along an edge
     */
    template < typename DeviceType, typename PotentialAccessorType >
    struct electric_flux_on_facet
    {
      typedef typename DeviceType::mesh_type MeshType;

      electric_flux_on_facet(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) {}

      std::vector<double> operator()(viennagrid_element_id facet) const
      {
        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim));

        viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device_.mesh(), facet, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

        std::vector<double> ret(3);
        if (cells_on_facet_begin + 1 == cells_on_facet_end)
          return ret;

        viennagrid_element_id c1 = cells_on_facet_begin[0];
        viennagrid_element_id c2 = cells_on_facet_begin[1];

        viennashe::permittivity_accessor<DeviceType>  permittivity(device_);

        const double potential_center = potential_.get_value(c1);
        const double potential_outer  = potential_.get_value(c2);

        std::vector<double> centroid_cell(3);
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c1, &(centroid_cell[0])));
        std::vector<double> centroid_other_cell(3);
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c2, &(centroid_other_cell[0])));
        double connection_len;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(3, &(centroid_cell[0]), &(centroid_other_cell[0]), &connection_len));

        std::vector<double> facet_centroid(3);
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), facet, &(facet_centroid[0])));
        double connection_in_cell;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(3, &(facet_centroid[0]), &(centroid_cell[0]), &connection_in_cell));
        double connection_in_other_cell;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(3, &(facet_centroid[0]), &(centroid_other_cell[0]), &connection_in_other_cell));

        const double permittivity_mean = (connection_in_cell + connection_in_other_cell) /
                                           (connection_in_cell/permittivity(c1) + connection_in_other_cell/permittivity(c2));

        const double Emag  = -(potential_outer - potential_center) / connection_len;

        ret[0] = Emag * permittivity_mean;
        return ret;
      }

      private:
        DeviceType            const & device_;
        PotentialAccessorType const & potential_;

    };

  } // namespace detail


  /**
   * @brief An electric flux accessor. Provides the electric flux along edges and on vertices given the electrostatic potential.
   */
  template < typename DeviceType, typename PotentialAccessorType >
  struct electric_flux_wrapper
  {
  public:
    typedef std::vector<double> value_type;

    electric_flux_wrapper(DeviceType const & device, PotentialAccessorType const & potential)
      : device_(device), potential_(potential)
    { }

    std::vector<double> operator()(viennagrid_element_id cell_or_facet) const
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim));

      viennashe::detail::electric_flux_on_facet<DeviceType, PotentialAccessorType> facet_eval(device_, potential_);

      std::vector<double> ret(3);
      if (viennagrid_topological_dimension_from_element_id(cell_or_facet) == cell_dim) // cell
      {
        viennashe::materials::checker no_conductor_filter(MATERIAL_NO_CONDUCTOR_ID);

        if (!no_conductor_filter(cell_or_facet))
          return ret;

        viennashe::util::value_holder_functor<std::vector<double> > result;

        viennashe::util::dual_box_flux_to_cell(device_,
                                               cell_or_facet,
                                               result, facet_eval);
        return result();

      }
      else if (viennagrid_topological_dimension_from_element_id(cell_or_facet) == cell_dim - 1) // facet
      {
        ret = facet_eval(cell_or_facet);
      }
      else
        throw std::runtime_error("current_density_wrapper::operator(): invalid element dimension!");

      return ret;


    }

    private:
      DeviceType            const & device_;
      PotentialAccessorType const & potential_;

  }; // electric_flux_wrapper


} // viennashe

#endif /* VIENNASHE_DD_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP */

