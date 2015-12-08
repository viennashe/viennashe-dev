#ifndef VIENNASHE_MAPPING_HPP
#define VIENNASHE_MAPPING_HPP

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
#include <fstream>
#include <vector>
#include <memory>

//ViennaGrid includes:
#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"
#include "viennashe/materials/all.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/she/timestep_quantities.hpp"
#include "viennashe/she/assemble_common.hpp"

#include "viennashe/util/misc.hpp"

/** @file viennashe/mapping.hpp
    @brief Distributes the unknown indices ('mapping') for the Poisson equation and the continuity equations
 */

namespace viennashe
{

  ////////////// Spatial quantity mapping ///////////////////////

  /**
   * @brief Maps unkown indices in the system matrix to the given unknown spatial quantity
   *
   * @param device The device on which the quantity is defined
   * @param quantity The unkown quantity
   * @param start_index The index offset, usefull if you have more than one unkown quantity
   * @return The first mapping index
   */
  template <typename VertexT>
  long create_mapping(viennashe::device const & device,
                      viennashe::unknown_quantity<VertexT> & quantity,
                      long start_index = 0)
  {
    //
    // Run the mapping on allowed vertices
    //
    long mapping_index = start_index;

    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      if ( quantity.get_boundary_type(*cit) != BOUNDARY_DIRICHLET && quantity.get_unknown_mask(*cit) )  // Dirichlet boundary condition
        quantity.set_unknown_index(*cit, mapping_index++);
      else
        quantity.set_unknown_index(*cit, -1);
    }

    return mapping_index;
  } //create_mapping


  ////////////// SHE quantity mapping ////////////////////

  namespace detail
  {

    /**
     * @brief Handles the mapping of a SHE unkown quantity on a single cell and H-space index
     * @param device   The device ...
     * @param cell     The Cell for which to do the mapping
     * @param index_H  The H index to map
     * @param quan     The quantity the mapping is for
     * @param conf     The simulator configuration
     * @param current_index The current index
     */
    template <typename CellType, typename SHEQuantity>
    void map_cell(viennashe::device const & device,
                  CellType const & cell,
                  std::size_t index_H,
                  SHEQuantity & quan,
                  viennashe::config const & conf,
                  long & current_index)
    {
      const long   expansion_order    = static_cast<long>(quan.get_expansion_order(cell, index_H));
      const double kinetic_energy_max = conf.max_kinetic_energy_range(quan.get_carrier_type_id());

      if ( device.get_material(cell) == viennashe::materials::metal::id) //this is a contact element attached to a semiconductor
      {
        if (viennashe::she::averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), cell, index_H) > 0
            && quan.get_kinetic_energy(cell, index_H) < kinetic_energy_max) // truncate range for high energies
        {
          quan.set_unknown_index(cell, index_H, current_index);
          current_index += viennashe::she::even_unknowns_on_node(expansion_order);
        }
        else
        {
          quan.set_unknown_index(cell, index_H, -1);
          quan.set_expansion_order(cell, index_H, 0);
        }
      }
      else if ( viennashe::she::averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), cell, index_H) > 0
                && quan.get_kinetic_energy(cell, index_H) < kinetic_energy_max) // truncate range for high energies
      {
        quan.set_unknown_index(cell, index_H, current_index);
        current_index += viennashe::she::even_unknowns_on_node(expansion_order);
      }
      else
      {
        quan.set_unknown_index(cell, index_H, -1);
        quan.set_expansion_order(cell, index_H, 0);
      }
    }


    /** @brief Assigns an unknown index for odd SHE unknowns to a facet.
     *  Preconditions:
     *   - Even unknowns on vertices already assigned
     */
    template <typename SHEQuantity>
    void map_facet(viennashe::device const & device,
                   SHEQuantity & quan,
                   viennashe::config const & conf,
                   viennagrid_element_id facet,
                   std::size_t index_H,
                   long & current_index)
    {
      typedef typename viennashe::device::mesh_type           MeshType;

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      // Get cells of the facet
      viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), facet, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

      if (cells_on_facet_begin + 1 == cells_on_facet_end)
      {
        quan.set_unknown_index(facet, index_H, -1);
        quan.set_expansion_order(facet, index_H, 0);
      }
      else
      {
        viennagrid_element_id c1 = cells_on_facet_begin[0];
        viennagrid_element_id c2 = cells_on_facet_begin[1];

        //
        // assign unknowns to facet if both cells carry unknowns:
        //
        long dof1 = quan.get_unknown_index(c1, index_H);
        long dof2 = quan.get_unknown_index(c2, index_H);
        if (dof1 > -1 && dof2 > -1
            && (averaged_density_of_states(quan, conf.dispersion_relation(quan.get_carrier_type_id()), facet, index_H) > 0)
            && (viennashe::materials::is_semiconductor(device.get_material(c1)) || viennashe::materials::is_semiconductor(device.get_material(c2))) ) // at least one of the two cells must be a semiconductor!
        {
          quan.set_unknown_index(facet, index_H, current_index);
          current_index += static_cast<long>(quan.get_unknown_num(facet, index_H));
        }
        else //energy is in forbidden band
        {
          quan.set_unknown_index(facet, index_H, -1);
          quan.set_expansion_order(facet, index_H, 0);
        }
      }

    } //map_edge()

  } //namespace detail




  /** @brief Distributes SHE unknown indices (dofs) over the device.
   *
   * @param device        The device (including a ViennaGrid mesh) on which simulation is carried out
   * @param quan          The quantity for which to create the mapping
   * @param conf          The simulator configuration
   * @param unknown_offset The first index to be assigned (allows for offsets); Nonzero value if e.g. the potential is also considered within Newton iteration
   */
  template <typename SHEQuantity>
  long create_even_mapping(viennashe::device const & device,
                           SHEQuantity & quan,
                           viennashe::config const & conf,
                           long unknown_offset = 0)
  {
    //
    // Distribute SHE unknown indices over the mesh.
    // This is NOT completely arbitrary:
    //  * First all even unknowns for a certain energy level
    //  * Then all even unknowns for the next energy level, and so on, until all energies for even are populated
    //  * Now the same for all even unknowns of other quantities when using Newton's method (call this function again)
    //  * Odd unknowns get mapped only after all even unknowns and spatial unknowns (potential, etc.) have been mapped.
    // Reasons:
    //  * after elimination of odd unknowns, "upper left block" of system matrix remains.
    //  * for each energy level, preconditioners for the "upper left block" can be computed in parallel
    //
    // Note that for a particular energy, only the first unknown-index (for Y_00) is stored.
    // Number of unknowns at that energy point is obtained from the (local) expansion order.
    //

    long unknown_index = unknown_offset;

    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));

    for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
    {
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (quan.get_unknown_mask(*cit) == false )
        {
          //no SHE unknowns here
            quan.set_unknown_index(*cit, index_H, -1);
        }
        else
        {
          if ( (index_H == 0) || (index_H == quan.get_value_H_size() - 1) )
            //no SHE unknowns at energy boundary (otherwise scattering at simulation domain gets pretty delicate w.r.t. symmetrization of in- and out-scattering)
            quan.set_unknown_index(*cit, index_H, -1);
          else
            detail::map_cell(device, *cit, index_H, quan, conf, unknown_index);
        }
      }
    }

    return unknown_index;
  }



  /** @brief Distributes SHE unknown indices (dofs) over the device.
   *
   * @param device        The device (including a ViennaGrid mesh) on which simulation is carried out
   * @param quan          The quantity for which to create the mapping
   * @param conf          The simulator configuration
   * @param unknown_offset The first index to be assigned (allows for offsets)
   */
  template <typename SHEQuantity>
  long create_odd_mapping(viennashe::device const & device,
                          SHEQuantity & quan,
                          viennashe::config const & conf,
                          long unknown_offset = 0)  //nonzero value if e.g. the potential is also considered within Newton iteration
  {
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *facets_begin, *facets_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim - 1, &facets_begin, &facets_end));


    long unknown_index = unknown_offset;

    for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
    {
      for (viennagrid_element_id *fit  = facets_begin;
                                  fit != facets_end;
                                ++fit)
      {
        detail::map_facet(device, quan, conf, *fit, index_H, unknown_index);
      }
    }

    return unknown_index;
  }




  /////////////// Public interface //////////////////

  /**
   * @brief Creates an unkown mapping for spatial quantities
   * @param device The device on which the quantities are defined
   * @param quantities A list of quantities
   * @param conf The simulator configuration
   * @return Information on the number of unknowns per quantity
   */
  inline map_info_type create_mapping(viennashe::device const & device,
                                      viennashe::she::timestep_quantities & quantities,
                                      viennashe::config const & conf)
  {
    typedef typename map_info_type::value_type::second_type   PairType;

    const bool use_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

    long last_mapping_index = 0;
    long current_mapping_index = 0;
    map_info_type map_info;

    // Create mapping for spatial quantities:
    for (std::size_t i=0; i<quantities.unknown_quantities().size(); ++i)
    {
      current_mapping_index  = viennashe::create_mapping(device, quantities.unknown_quantities()[i], last_mapping_index);

      map_info[quantities.unknown_quantities()[i].get_name()] = PairType(std::size_t(current_mapping_index - last_mapping_index), 0);
      if (use_newton)
        last_mapping_index = current_mapping_index;
    }

    // Create mapping for SHE quantities (even):
    for (std::size_t i=0; i<quantities.unknown_she_quantities().size(); ++i)
    {
      current_mapping_index  = viennashe::create_even_mapping(device, quantities.unknown_she_quantities()[i], conf, last_mapping_index);

      map_info[quantities.unknown_she_quantities()[i].get_name()].first = std::size_t(current_mapping_index - last_mapping_index);

      if (use_newton)
        last_mapping_index = current_mapping_index;
      else //mapping for odd quantities
      {
        long last_even_mapping_index = current_mapping_index;
        current_mapping_index  = viennashe::create_odd_mapping(device, quantities.unknown_she_quantities()[i], conf, last_even_mapping_index);

        map_info[quantities.unknown_she_quantities()[i].get_name()].second = std::size_t(current_mapping_index - last_even_mapping_index);
      }
    }

    // Create mapping for SHE quantities (odd) when using Newton:
    if (use_newton)
    {
      for (std::size_t i=0; i<quantities.unknown_she_quantities().size(); ++i)
      {
        current_mapping_index  = viennashe::create_odd_mapping(device, quantities.unknown_she_quantities()[i], conf, last_mapping_index);

        map_info[quantities.unknown_she_quantities()[i].get_name()].second = std::size_t(current_mapping_index - last_mapping_index);
        last_mapping_index = current_mapping_index;
      }
    }

    return map_info;
  } // create_mapping()

} //namespace viennashe

#endif
