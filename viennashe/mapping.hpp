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
#include <petsc.h>

//ViennaGrid includes:
#include "viennagrid/mesh/mesh.hpp"

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
  template <typename DeviceType, typename VertexT>
  long create_mapping(DeviceType const & device,
                      viennashe::unknown_quantity<VertexT> & quantity,
                      long start_index = 0)
  {
    typedef typename DeviceType::mesh_type           MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    //
    // Run the mapping on allowed vertices
    //
    long mapping_index = start_index;

    CellContainer cells(device.mesh());
    for (CellIterator cit = cells.begin();
        cit != cells.end();
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
    template <typename DeviceType, typename CellType, typename SHEQuantity>
    void map_cell(DeviceType const & device,
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
    template <typename DeviceType, typename SHEQuantity, typename FacetType>
    void map_facet(DeviceType const & device,
                   SHEQuantity & quan,
                   viennashe::config const & conf,
                   FacetType const & facet,
                   std::size_t index_H,
                   long & current_index)
    {
      typedef typename DeviceType::mesh_type           MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type            CellType;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

      // Get cells of the facet
      CellOnFacetContainer cells_on_facet(device.mesh(), viennagrid::handle(device.mesh(), facet));
      CellOnFacetIterator cofit = cells_on_facet.begin();

      CellType const & c1 = *cofit;
      ++cofit;
      if (cofit == cells_on_facet.end())
      {
        quan.set_unknown_index(facet, index_H, -1);
        quan.set_expansion_order(facet, index_H, 0);
      }
      else
      {
        CellType const & c2 = *cofit;

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
  template <typename DeviceType, typename SHEQuantity>
  long create_even_mapping(DeviceType const & device,
                           SHEQuantity & quan,
                           viennashe::config const & conf,
                           long unknown_offset = 0)
  {
    typedef typename DeviceType::mesh_type              MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    MeshType const & mesh = device.mesh();

    CellContainer cells(mesh);

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

    /*
        * PETSC
        * */
    if(quan.get_value_H_size() == 0 )return 0;

    int size = 1, rank = 0;
    if (conf.linear_solver().id() == viennashe::solvers::linear_solver_config::petsc_parallel_linear_solver)
    {
      MPI_Comm_size(PETSC_COMM_WORLD,&size);
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    }
    std::size_t maxquantityfull = quan.get_value_H_size();
    //* Overleap the regions to avoid the singular matrix*/
    std::size_t maxquantity = rank < size-1 ? ((rank+1)*maxquantityfull)/size+1 :maxquantityfull ;
    std::size_t init = ((rank)*maxquantityfull)/size;//;
    long unknown_index = unknown_offset;
    for (std::size_t index_H = init; index_H < maxquantity; ++index_H)
    {
      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
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
  template <typename DeviceType, typename SHEQuantity>
  long create_odd_mapping(DeviceType const & device,
                          SHEQuantity & quan,
                          viennashe::config const & conf,
                          long unknown_offset = 0)  //nonzero value if e.g. the potential is also considered within Newton iteration
  {
    typedef typename DeviceType::mesh_type              MeshType;

    typedef typename viennagrid::result_of::const_facet_range<MeshType>::type         FacetContainer;
    typedef typename viennagrid::result_of::iterator<FacetContainer>::type            FacetIterator;

    MeshType const & mesh = device.mesh();

    FacetContainer facets(mesh);
    /*
     * PETSC
     * */
    if(quan.get_value_H_size() == 0 )return 0;
    int size = 1,rank = 0;
    if (conf.linear_solver().id() == viennashe::solvers::linear_solver_config::petsc_parallel_linear_solver)
    {
      MPI_Comm_size(PETSC_COMM_WORLD,&size);
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    }
    std::size_t maxquantityfull = quan.get_value_H_size();
    std::size_t maxquantity = rank < size-1 ? ((rank+1)*maxquantityfull)/size+1:maxquantityfull ;
    std::size_t init = ((rank)*maxquantityfull)/size;
    long unknown_index = unknown_offset;
    for (std::size_t index_H = init; index_H < maxquantity; ++index_H)
    {
      for (FacetIterator fit = facets.begin();
                         fit != facets.end();
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
  template<typename DeviceT>
  map_info_type create_mapping(DeviceT const & device,
                               viennashe::she::timestep_quantities<DeviceT> & quantities,
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
