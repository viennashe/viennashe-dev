#ifndef VIENNASHE_SHE_BOUNDARY_CONDITIONS_HPP
#define VIENNASHE_SHE_BOUNDARY_CONDITIONS_HPP

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
#include <limits>

// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/util/misc.hpp"
#include "viennashe/forwards.h"
#include "viennashe/accessors.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/config.hpp"
#include "viennashe/she/she_quantity.hpp"
#include "viennashe/she/timestep_quantities.hpp"

/** @file viennashe/she/boundary_conditions.hpp
    @brief Writes the SHE boundary conditions to the mesh
*/


namespace viennashe
{
  namespace she
  {

    /** @brief Writes boundary conditions for SHE to the device. Stores the result using ViennaData.
     *
     * @param device  The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan    The unknown SHE quantities on edges and vertices
     * @param conf    The simulator configuration
     */
    template <typename DeviceType,
              typename VertexT, typename EdgeT>
    void write_boundary_conditions(DeviceType const & device,
                                   viennashe::she::unknown_she_quantity<VertexT, EdgeT> & quan,
                                   viennashe::config const & conf)
    {
      viennagrid_mesh mesh = device.mesh();

      viennashe::contact_carrier_density_accessor<DeviceType> bnd_carrier_density(device, quan.get_carrier_type_id());

      const double kB = viennashe::physics::constants::kB;
      viennashe::config::dispersion_relation_type dispersion = conf.dispersion_relation(quan.get_carrier_type_id());

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        const double T = device.get_lattice_temperature(*cit);

        // Step 1: Find normalization
        //
        // Compute boundary correction such that boundary values are given by: n_conc * exp(-eps/(kB*T)) / bnd_corr
        //
        double bnd_corr = 0;

        for (std::size_t index_H = 1;
                         index_H < quan.get_value_H_size() - 1;
                       ++index_H)
        {
          double energy_mid = quan.get_kinetic_energy(*cit, index_H);
          double energy_lower = std::max<double>((quan.get_kinetic_energy(*cit, index_H - 1) + energy_mid) / 2.0, 0.0);
          double energy_upper = (quan.get_kinetic_energy(*cit, index_H + 1) + energy_mid) / 2.0;

          if (energy_lower < energy_upper && (energy_upper >= 0))
          {
            double height = box_height(quan, *cit, index_H);

            bnd_corr += exp(-energy_mid / (kB*T)) * averaged_density_of_states(quan, dispersion, *cit, index_H) * height;
          }
        }
        bnd_corr *= dispersion.symmetry_factor();

        //
        // Step 2: Write boundary values
        //
        for (std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          double energy = quan.get_kinetic_energy(*cit, index_H);
          double energy_lower = std::max<double>((quan.get_kinetic_energy(*cit, index_H - 1) + energy) / 2.0, 0.0);
          double energy_upper = (quan.get_kinetic_energy(*cit, index_H + 1) + energy) / 2.0;

          double bnd_value = 0.0;
          if (energy_lower < energy_upper && (energy_upper >= 0)) // ensure conservation of density upon integration
          {
            const double TL = device.get_lattice_temperature(*cit);

            switch (conf.she_discretization_type())
            {
            case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
              bnd_value = bnd_carrier_density(*cit, TL) * exp(-energy / (kB * T)) / bnd_corr;
              break;
            case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
              bnd_value = bnd_carrier_density(*cit, TL) * averaged_density_of_states(quan, dispersion, *cit, index_H) * exp(-energy / (kB * T)) / bnd_corr;
              break;
            default:
              throw std::runtime_error("Unknown SHE discretization type");
            }
          }
          quan.set_boundary_value(*cit, index_H, bnd_value);
        } // for index_H

      } // for cells

    } // write_boundary_conditions())


    template <typename DeviceType>
    void write_boundary_conditions(DeviceType const & device,
                                   timestep_quantities<DeviceType> & quantities,
                                   viennashe::config const & conf)
    {
      typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_she_quantity_type  SHEUnknownType;

      SHEUnknownType & f_n = quantities.electron_distribution_function();
      SHEUnknownType & f_p = quantities.hole_distribution_function();

      // Setup energies for electrons:
      if (conf.with_electrons() && conf.get_electron_equation() == EQUATION_SHE)
        write_boundary_conditions(device, f_n, conf);

      // Setup energies for holes:
      if (conf.with_holes() && conf.get_hole_equation() == EQUATION_SHE)
        write_boundary_conditions(device, f_p, conf);
    }

  }
}

#endif
