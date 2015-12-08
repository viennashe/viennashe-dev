#ifndef VIENNASHE_SHE_SETUP_ENERGIES_HPP
#define VIENNASHE_SHE_SETUP_ENERGIES_HPP

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
#include <limits>

#include "viennagrid/viennagrid.h"

#include "viennashe/util/misc.hpp"
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/config.hpp"

#include "viennashe/she/she_quantity.hpp"
#include "viennashe/she/timestep_quantities.hpp"
#include "viennashe/simulator_quantity.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/accessors.hpp"

/** @file viennashe/she/setup_energies.hpp
    @brief Writes the kinetic energies to each nodes in the (x,H)-domain
*/


namespace viennashe
{
  namespace she
  {


    /** @brief Computes the kinetic energy over the device as obtained from the provided potential.
     *
     * @param device        The device on which simulation is carried out
     * @param quan          The SHE quantities
     * @param conf          The simulator configuration
     * @param potential     The potential as computed with Poisson equation
     * @param quantum_correction The quantum correction potential
     */
    template <typename DeviceT>
    void setup_energies(DeviceT const & device,
                        viennashe::she::unknown_she_quantity<double> & quan,
                        viennashe::config const & conf,
                        viennashe::unknown_quantity<double> const & potential,
                        viennashe::unknown_quantity<double> const & quantum_correction)
    {
      viennagrid_mesh mesh = device.mesh();

      double max_potential = -1.0 * std::numeric_limits<double>::max();
      double min_potential =        std::numeric_limits<double>::max();

      //
      // Step 1: Find minimum of band edge center and write band edge center to device.
      //         Only consider vertices which are attached to semiconductor cells.
      //

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (!viennashe::materials::is_semiconductor(device.get_material(*cit)))
          continue;

        const double cell_potential = potential.get_value(*cit) + quantum_correction.get_value(*cit);

        if (cell_potential < min_potential)  min_potential = cell_potential;
        if (cell_potential > max_potential)  max_potential = cell_potential;
      }

      //
      // Step 2: Set total energy
      //

      double total_energy_range = viennashe::physics::constants::q * (max_potential - min_potential) + 2.0 * conf.energy_spacing() + conf.min_kinetic_energy_range(quan.get_carrier_type_id());
      double total_energy_bandedge = (quan.get_carrier_type_id() == ELECTRON_TYPE_ID) ? -viennashe::physics::constants::q * max_potential - conf.energy_spacing() + viennashe::materials::si::band_gap() / 2.0
                                                                                      : -viennashe::physics::constants::q * min_potential + conf.energy_spacing() - viennashe::materials::si::band_gap() / 2.0;
      double total_energy_increment = ((quan.get_carrier_type_id() == ELECTRON_TYPE_ID) ? 1.0 : -1.0) * conf.energy_spacing();

      if (!conf.use_h_transformation())
      {
        total_energy_range = conf.min_kinetic_energy_range(quan.get_carrier_type_id());
        total_energy_bandedge = -total_energy_increment; // align band edge with index_H=1
      }

      // Set total energy
      viennagrid_int facet_count;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(mesh, cell_dim - 1, &facet_count));
      std::size_t num_energies = static_cast<std::size_t>(total_energy_range / conf.energy_spacing()) + 1;
      quan.resize(cells_end - cells_begin, facet_count, num_energies);
      for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
      {
        quan.set_value_H(index_H, total_energy_bandedge + index_H * total_energy_increment);
      }

      // Set bandedge shift on vertices and edges (only when using H-transform):
      if (conf.use_h_transformation())
      {
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          if (quan.get_carrier_type_id() == ELECTRON_TYPE_ID)
            quan.set_bandedge_shift(*cit, -viennashe::physics::constants::q * potential.get_value(*cit) + viennashe::materials::si::band_gap() / 2.0);
          else
            quan.set_bandedge_shift(*cit, -viennashe::physics::constants::q * potential.get_value(*cit) - viennashe::materials::si::band_gap() / 2.0);
        }

        viennagrid_element_id *facets_begin, *facets_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim - 1, &facets_begin, &facets_end));
        for (viennagrid_element_id *fit = facets_begin;
                                    fit != facets_end;
                                  ++fit)
        {
          viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

          if (cells_on_facet_begin + 1 == cells_on_facet_end)
            quan.set_bandedge_shift(*fit, quan.get_bandedge_shift(cells_on_facet_begin[0]));
          else
            quan.set_bandedge_shift(*fit, (quan.get_bandedge_shift(cells_on_facet_begin[0]) + quan.get_bandedge_shift(cells_on_facet_begin[1])) / 2.0);
        }
      }

    }

    /** @brief Computes the kinetic energy over the device as obtained from the provided potential.
     *
     * @param device        The device on which simulation is carried out
     * @param quantities    The quantities defined on the device
     * @param conf          The simulator configuration
     */
    inline void setup_energies(viennashe::device const & device,
                               timestep_quantities & quantities,
                               viennashe::config const & conf)
    {
      typedef typename viennashe::she::timestep_quantities::unknown_quantity_type      SpatialUnknownType;
      typedef typename viennashe::she::timestep_quantities::unknown_she_quantity_type  SHEUnknownType;

      SpatialUnknownType const & potential      = quantities.get_unknown_quantity(viennashe::quantity::potential());
      SpatialUnknownType const & quantum_corr_n = quantities.get_unknown_quantity(viennashe::quantity::density_gradient_electron_correction());
      SpatialUnknownType const & quantum_corr_p = quantities.get_unknown_quantity(viennashe::quantity::density_gradient_hole_correction());
      SHEUnknownType           & f_n            = quantities.electron_distribution_function();
      SHEUnknownType           & f_p            = quantities.hole_distribution_function();

      // Setup energies for electrons:
      if (conf.with_electrons() && conf.get_electron_equation() == EQUATION_SHE)
        setup_energies(device, f_n, conf, potential, quantum_corr_n);

      // Setup energies for holes:
      if (conf.with_holes() && conf.get_hole_equation() == EQUATION_SHE)
        setup_energies(device, f_p, conf, potential, quantum_corr_p);
    }

  }
}

#endif
