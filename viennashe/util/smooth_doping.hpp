#ifndef VIENNASHE_UTIL_SMOOTH_DOPING_HPP
#define VIENNASHE_UTIL_SMOOTH_DOPING_HPP

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

/** @file viennashe/util/smooth_doping.hpp
    @brief Routines for ensuring a constant doping along contacts.
*/

#include <set>
#include <queue>
#include <vector>

#include "viennagrid/viennagrid.h"

#include "viennashe/device.hpp"
#include "viennashe/materials/all.hpp"

namespace viennashe
{
  namespace detail
  {
    /** @brief Utility function for smoothing the doping next to contacts.
      *
      * Otherwise, all types of funny (and not so funny) effects can happen, potentially making the simulation result complete bogus.
      */
    template <typename DeviceT>
    void smooth_doping_at_contacts(DeviceT & d)
    {
      viennagrid_mesh mesh = d.mesh();

      viennagrid_dimension cell_dim;
      viennagrid_mesh_cell_dimension_get(mesh, &cell_dim);

      viennagrid_element_id *cells_begin, *cells_end;
      viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end);

      //
      // Step 1: Enumerate contacts (not relying on segments)
      //

      //
      // Enumerates the contacts starting from value 1.
      // The first contact cell encountered is assigned a value of 1, which is propagated to all its neighbors (reachable via other conductor cells).
      // Second contact cell not part of first contact is assigned 2, followed by propagation to all other cells of the same contact.
      // Same for other contacts.
      // Non-contact cells remain at a contact ID 0.
      //
      std::vector<std::size_t> contact_id(cells_end - cells_begin);
      std::set<viennagrid_element_id> cells_semiconductor;

      std::size_t current_contact_id = 1;

      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (viennashe::materials::is_conductor(d.get_material(*cit)))
        {
          std::size_t cell_index = std::size_t(viennagrid_index_from_element_id(*cit));

          // Only consider further if contact ID not yet assigned
          if (contact_id[cell_index] == 0)
          {
            contact_id[cell_index] = current_contact_id;

            std::queue<viennagrid_element_id> cells_queue;
            std::set<viennagrid_element_id> cells_discovered;

            cells_queue.push(*cit);
            cells_discovered.insert(*cit);

            //
            // Iteration over all conductor cells attached to this cell (breadth-first iteration)
            //
            while (!cells_queue.empty())
            {
              viennagrid_element_id current_cell = cells_queue.front();
              cells_queue.pop();

              viennagrid_element_id *neighbors_begin, *neighbors_end;
              viennagrid_element_neighbor_elements(mesh, current_cell, cell_dim - 1, cell_dim, &neighbors_begin, &neighbors_end);
              for (viennagrid_element_id *nit = neighbors_begin; nit != neighbors_end; ++nit)
              {
                // conductor cells are added to the queue to be revisited later
                if (viennashe::materials::is_conductor(d.get_material(*nit)))
                {
                  if (cells_discovered.find(*nit) == cells_discovered.end())
                  {
                    cells_discovered.insert(*nit);
                    cells_queue.push(*nit);
                    contact_id[viennagrid_index_from_element_id(*nit)] = current_contact_id;
                  }
                }
                else if (viennashe::materials::is_semiconductor(d.get_material(*nit)))
                {
                  // Remember this cell for later averaging of the doping
                  cells_semiconductor.insert(*nit);

                  if (contact_id[viennagrid_index_from_element_id(*nit)] != 0 && contact_id[viennagrid_index_from_element_id(*nit)] != current_contact_id)
                    log::warning() << "Warning: Cell " << viennagrid_index_from_element_id(*nit) << " is attached to at least two different contacts. Doping smoother might give inconsistent results!" << std::endl;

                  contact_id[viennagrid_index_from_element_id(*nit)] = current_contact_id;
                }

              } //for
            } // while

            // all (multi-hop) conducting neighbors tagged, increase contact_id:
            ++current_contact_id;

          } // if no contact ID assigned yet
        } // if conductor
      } // for cells


      //
      // Step 2: Compute doping averages over semiconductor cells attached to contact cells:
      //
      std::vector<double> avg_doping_n_values(current_contact_id-1);
      std::vector<double> avg_doping_p_values(current_contact_id-1);

      for (typename std::set<viennagrid_element_id>::const_iterator it  = cells_semiconductor.begin();
                                                                    it != cells_semiconductor.end();
                                                                  ++it)
      {
        std::size_t corrected_contact_id = contact_id[viennagrid_index_from_element_id(*it)] - 1;
        avg_doping_n_values[corrected_contact_id] = std::max(d.get_doping_n(*it), avg_doping_n_values[corrected_contact_id]);
        avg_doping_p_values[corrected_contact_id] = std::max(d.get_doping_p(*it), avg_doping_p_values[corrected_contact_id]);
      }


      //
      // Step 3: Set new doping:
      //
      for (typename std::set<viennagrid_element_id>::const_iterator it  = cells_semiconductor.begin();
                                                                    it != cells_semiconductor.end();
                                                                  ++it)
      {
        std::size_t corrected_contact_id = contact_id[viennagrid_index_from_element_id(*it)] - 1;
        d.set_doping_n(avg_doping_n_values[corrected_contact_id], *it);
        d.set_doping_p(avg_doping_p_values[corrected_contact_id], *it);
      }

    }
  }


} //namespace viennashe


#endif
