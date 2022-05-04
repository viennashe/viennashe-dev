#ifndef VIENNASHE_UTIL_SMOOTH_DOPING_HPP
#define VIENNASHE_UTIL_SMOOTH_DOPING_HPP

/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
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

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/mesh/neighbor_iteration.hpp"

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
      typedef typename DeviceT::mesh_type   MeshType;

      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

      typedef typename viennagrid::result_of::const_neighbor_range<MeshType, CellType, FacetType>::type    NeighborRange;
      typedef typename viennagrid::result_of::iterator<NeighborRange>::type                                NeighborIterator;

      MeshType const & mesh = d.mesh();

      CellContainer cells(mesh);

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
      std::vector<std::size_t> contact_id(cells.size());
      std::set<std::size_t> cells_semiconductor;

      std::size_t current_contact_id = 1;

      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        if (viennashe::materials::is_conductor(d.get_material(*cit)))
        {
          // Only consider further if contact ID not yet assigned
          if (contact_id[std::size_t(cit->id().get())] == 0)
          {
            contact_id[std::size_t(cit->id().get())] = current_contact_id;

            std::queue<std::size_t> cells_queue;
            std::set<std::size_t> cells_discovered;

            cells_queue.push(std::size_t(cit->id().get()));
            cells_discovered.insert(std::size_t(cit->id().get()));

            //
            // Iteration over all conductor cells attached to this cell (breadth-first iteration)
            //
            while (!cells_queue.empty())
            {
              std::size_t current_cell_id = cells_queue.front();
              cells_queue.pop();

              CellType const & current_cell = cells[current_cell_id];

              NeighborRange neighbors(mesh, viennagrid::handle(mesh, current_cell));
              for (NeighborIterator nit = neighbors.begin(); nit != neighbors.end(); ++nit)
              {
                std::size_t neighbor_id = std::size_t(nit->id().get());

                // conductor cells are added to the queue to be revisited later
                if (viennashe::materials::is_conductor(d.get_material(*nit)))
                {
                  if (cells_discovered.find(neighbor_id) == cells_discovered.end())
                  {
                    cells_discovered.insert(neighbor_id);
                    cells_queue.push(neighbor_id);
                    contact_id[neighbor_id] = current_contact_id;
                  }
                }
                else if (viennashe::materials::is_semiconductor(d.get_material(*nit)))
                {
                  // Remember this cell for later averaging of the doping
                  cells_semiconductor.insert(neighbor_id);

                  if (contact_id[neighbor_id] != 0 && contact_id[neighbor_id] != current_contact_id)
                    log::warning() << "Warning: Cell " << neighbor_id << " is attached to at least two different contacts. Doping smoother might give inconsistent results!" << std::endl;

                  contact_id[neighbor_id] = current_contact_id;
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

      for (typename std::set<std::size_t>::const_iterator it  = cells_semiconductor.begin();
                                                          it != cells_semiconductor.end();
                                                        ++it)
      {
        std::size_t cell_id = *it;
        CellType const & cell = cells[cell_id];

        std::size_t corrected_contact_id = contact_id[cell_id] - 1;
        avg_doping_n_values[corrected_contact_id] = std::max(d.get_doping_n(cell), avg_doping_n_values[corrected_contact_id]);
        avg_doping_p_values[corrected_contact_id] = std::max(d.get_doping_p(cell), avg_doping_p_values[corrected_contact_id]);
      }


      //
      // Step 3: Set new doping:
      //
      for (typename std::set<std::size_t>::const_iterator it  = cells_semiconductor.begin();
                                                          it != cells_semiconductor.end();
                                                        ++it)
      {
        std::size_t cell_id = *it;
        CellType const & cell = cells[cell_id];

        std::size_t corrected_contact_id = contact_id[cell_id] - 1;
        d.set_doping_n(avg_doping_n_values[corrected_contact_id], cell);
        d.set_doping_p(avg_doping_p_values[corrected_contact_id], cell);
      }

    }
  }


} //namespace viennashe


#endif
