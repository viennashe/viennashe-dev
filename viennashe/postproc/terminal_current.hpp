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
                              SegmentT const & semiconductor,
                              SegmentT const & terminal)
  {
    double current = 0;

    /*
    CellContainer cells(terminal);
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      PointType centroid_cell = viennagrid::centroid(*cit);

      FacetOnCellContainer facets_on_cell(*cit);
      for (FacetOnCellIterator focit = facets_on_cell.begin();
          focit != facets_on_cell.end();
          ++focit)
      {
        if ( viennagrid::is_interface(terminal, semiconductor, *focit) )
        {
          CellOnFacetContainer cells_on_facet(device.mesh(), focit.handle());

          CellType const *other_cell_ptr = util::get_other_cell_of_facet(device.mesh(), *focit, *cit);

          if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

          PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
          PointType cell_connection = centroid_other_cell - centroid_cell;
          PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
          PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);
          const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);

          if ( &(*cit) == &(cells_on_facet[0]) )
            current += current_accessor(*focit) * weighted_interface_area;
          else  //reference direction is opposite of what we need
            current -= current_accessor(*focit) * weighted_interface_area;

          //log::info() << *eovit << std::endl;
          //log::info() << current_density_accessor( *eovit ) << " * " << effective_interface << " = " << current << std::endl;
        }// for edges
      }
    } */ // for vertices

    throw std::runtime_error("get_terminal_current(): TODO: Port to ViennaGrid 3.0");

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
  template<typename DeviceT, typename T1, typename T2, typename SegmentT>
  double get_terminal_current(DeviceT const & device,
                              viennashe::config const & conf,
                              viennashe::she::unknown_she_quantity<T1, T2> const & quan,
                              SegmentT const & semiconductor,
                              SegmentT const & terminal
                              )
  {
    typedef viennashe::she::unknown_she_quantity<T1, T2>  SHEQuantityType;

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

