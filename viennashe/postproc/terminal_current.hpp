#ifndef VIENNASHE_UTIL_TERMINAL_CURRENT_HPP
#define VIENNASHE_UTIL_TERMINAL_CURRENT_HPP

/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:      MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include "viennagrid/viennagrid.h"

#include "viennashe/util/filter.hpp"
#include "viennashe/accessors.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/she/she_quantity.hpp"
#include "viennashe/config.hpp"
#include "viennashe/she/postproc/current_density.hpp"
#include "viennashe/postproc/current_density.hpp"

/** @file viennashe/postproc/terminal_current.hpp
    @brief Convenience functions to access the terminal (contact) current
*/


namespace viennashe
{

  /**
   * @brief Returns the terminal current for a number of given vertices. Considers carrier flux and displacement currents.
   * @param device            The device
   * @param current_accessor  An accessor for the carrier current density
   * @param terminal          The terminal segment from which the current should be extracted
   * @param semiconductor     The semiconductor segment to which the current is flowing
   * @return The current in Ampere
   */
  template<typename DeviceT, typename CurrentAccessorT, typename SegmentT>
  double get_terminal_current(DeviceT const & device,
                              CurrentAccessorT const & current_accessor,
                              viennagrid_region semiconductor,
                              viennagrid_region terminal)
  {
    double current = 0;

    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      viennagrid_bool cell_in_terminal;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(terminal, *cit, &cell_in_terminal));

      if (!cell_in_terminal)
        continue;

      std::vector<viennagrid_numeric> centroid_cell(3);
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device.mesh(), *cit, &(centroid_cell[0])));

      viennagrid_element_id *facets_on_cell_begin, *facets_on_cell_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(device.mesh(), *cit, cell_dim - 1, &facets_on_cell_begin, &facets_on_cell_end));

      for (viennagrid_element_id *focit  = facets_on_cell_begin;
                                  focit != facets_on_cell_end;
                                ++focit)
      {
        viennagrid_bool facet_in_semiconductor;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(semiconductor, *cit, &facet_in_semiconductor));

        if (facet_in_semiconductor) // facet at semiconductor-terminal interface
        {
          viennagrid_element_id *other_cell_ptr;
          util::get_other_cell_of_facet(device.mesh(), *focit, *cit, &other_cell_ptr);

          if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

          std::vector<double> centroid_other_cell(3);
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device.mesh(), *other_cell_ptr, &(centroid_other_cell[0])));
          std::vector<double> cell_connection(3);
          for (std::size_t i=0; i<cell_connection.size(); ++i)
            cell_connection[i] = centroid_other_cell[i] - centroid_cell[i];
          double connection_len;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_norm_2(3, &(cell_connection[0]), &connection_len));
          std::vector<double> cell_connection_normalized(3);
          for (std::size_t i=0; i<cell_connection.size(); ++i)
            cell_connection_normalized[i] = cell_connection[i] / connection_len;
          std::vector<double> facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(device.mesh(), *cit, *focit);

          double facet_area;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_volume(device.mesh(), *focit, &facet_area));
          double facet_contribution;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_inner_prod(3, &(facet_unit_normal[0]), &(cell_connection_normalized[0]), &facet_contribution));
          const double weighted_interface_area = facet_area * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);

          viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), *focit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

          if (*cit == cells_on_facet_begin[0]) //reference direction is into box
            current += current_accessor(*focit) * weighted_interface_area;
          else  //reference direction is opposite of what we need
            current -= current_accessor(*focit) * weighted_interface_area;

          //log::info() << *eovit << std::endl;
          //log::info() << current_density_accessor( *eovit ) << " * " << effective_interface << " = " << current << std::endl;
        }// for edges
      }
    }  // for cells

    return current;
  } // get_terminal_current

  /**
   * @brief Returns the the terminal current for a given contact segment. Considers carrier flux and displacement currents.
   * @param device         The device
   * @param conf           The simulator configuration
   * @param quan           The quantity from which the current should be extracted (here: electron or hole distribution function)
   * @param terminal       The terminal segment from which the current should be extracted
   * @param semiconductor  The semiconductor segment to which the current is flowing
   * @return The current in Ampere
   */
  template<typename DeviceT, typename SegmentT>
  double get_terminal_current(DeviceT const & device,
                              viennashe::config const & conf,
                              viennashe::she::unknown_she_quantity<double> const & quan,
                              SegmentT const & semiconductor,
                              SegmentT const & terminal
                              )
  {
    typedef viennashe::she::unknown_she_quantity<double>  SHEQuantityType;

    viennashe::she::current_density_wrapper<DeviceT, SHEQuantityType> current_wrapper(device, conf, quan);
    return get_terminal_current(device, current_wrapper, semiconductor, terminal);
  }


  /**
   * @brief Returns the current (drift diffusion) from a semiconductor into a conductor segment
   * @param device The device
   * @param ctype Electron current or hole current
   * @param potential An accessor to the elec. potential
   * @param carrier An accessor to the carrier density (must match with ctype!)
   * @param mobility_model The mobility model (just an accessor to the mobility)
   * @param semi The semiconductor segment
   * @param conductor The conductor segment adjacent to the semiconductor
   * @return The current in Ampere
   */
  template < typename DeviceType,
             typename PotentialAccessor,
             typename AccessorTypeCarrier,
             typename MobilityModel,
             typename SegmentType >
  double get_terminal_current(DeviceType const & device,
                              viennashe::carrier_type_id ctype,
                              PotentialAccessor const & potential,
                              AccessorTypeCarrier const & carrier,
                              MobilityModel const & mobility_model,
                              SegmentType const & semi,
                              SegmentType const & conductor)
  {
    viennashe::current_density_wrapper<DeviceType, PotentialAccessor, AccessorTypeCarrier, MobilityModel> Jfield(device, ctype, potential, carrier, mobility_model);
    return viennashe::get_terminal_current(device, Jfield, semi, conductor);
  }

} // viennashe


#endif

