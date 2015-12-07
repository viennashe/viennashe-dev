#ifndef VIENNASHE_DD_POSTPROC_CURRENT_HPP
#define VIENNASHE_DD_POSTPROC_CURRENT_HPP

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

#include "viennashe/accessors.hpp"
#include "viennashe/scharfetter_gummel.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/checks.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/util/misc.hpp"
#include "viennashe/util/dual_box_flux.hpp"

#include "viennashe/she/postproc/macroscopic.hpp"
#include "viennashe/simulator_quantity.hpp"

/** @file viennashe/postproc/current_density.hpp
    @brief Computes the current density (both Jn and Jp) after using a drift diffusion solution
 */

namespace viennashe
{

  namespace detail
  {
    /** @brief An accessor to the current density (drift diffusion only!) on edges */
    template <typename DeviceType,
              typename PotentialAccessorType,
              typename AccessorTypeCarrier,
              typename MobilityModel>
    struct current_density_on_facet
    {
      typedef std::vector<double> value_type;

      current_density_on_facet(DeviceType const & device,
                               viennashe::carrier_type_id ctype,
                               PotentialAccessorType const & potential,
                               AccessorTypeCarrier const & carrier,
                               MobilityModel const & mobility_model)
      : device_(device), carrier_type_id_(ctype), potential_(potential), carrier_(carrier), mobility_(mobility_model) { }

      value_type operator()(viennagrid_element_id facet) const
      {
        typename viennashe::contact_carrier_density_accessor<DeviceType> bnd_carrier_density(device_, carrier_type_id_);

        scharfetter_gummel flux_approximator(carrier_type_id_);

        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim));

        viennagrid_element_id *cells_begin, *cells_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device_.mesh(), facet, cell_dim, &cells_begin, &cells_end));

        value_type ret(3);

        if (cells_end == cells_begin + 1) // facet on surface
          return ret;

        viennagrid_element_id c1 = cells_begin[0];
        viennagrid_element_id c2 = cells_begin[1];

        if (viennashe::materials::is_insulator(device_.get_material(c1)) || viennashe::materials::is_insulator(device_.get_material(c2))) // no current into insulator
          return ret;

        // at least one of the two cells must be a semiconductor:
        if (!viennashe::materials::is_semiconductor(device_.get_material(c1)) && !viennashe::materials::is_semiconductor(device_.get_material(c2)))
          return ret;

        const double potential_center = potential_(c1);
        const double carrier_center   = viennashe::materials::is_conductor(device_.get_material(c1))
                                        ? bnd_carrier_density(c2, device_.get_lattice_temperature(c2))
                                        : carrier_.get_value(c1);
        const double T  = device_.get_lattice_temperature(c1);

        const double mobility = mobility_(c1, c2, potential_);

        std::vector<viennagrid_numeric> centroid_1(3);
        std::vector<viennagrid_numeric> centroid_2(3);

        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c1, &(centroid_1[0])));
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c2, &(centroid_2[0])));

        viennagrid_numeric distance;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(cell_dim, &(centroid_1[0]), &(centroid_2[0]), &distance));

        const double potential_outer = potential_(c2);
        const double carrier_outer   = viennashe::materials::is_conductor(device_.get_material(c2))
                                       ? bnd_carrier_density(c1, device_.get_lattice_temperature(c1))
                                       : carrier_.get_value(c2);

        if ( carrier_outer <= 0 || carrier_center <= 0 ) return ret;

        const double polarity = (carrier_type_id_ == viennashe::ELECTRON_TYPE_ID) ? -1.0 : 1.0;
        const double charge_flux = flux_approximator(carrier_center, carrier_outer, potential_center, potential_outer, distance, mobility, T);
        const double Jmag = polarity * mobility * charge_flux;

        ret[0] = Jmag;
        return ret;
      } // operator()

    private:
      DeviceType                 const & device_;
      viennashe::carrier_type_id         carrier_type_id_;
      PotentialAccessorType      const & potential_;
      AccessorTypeCarrier        const & carrier_;
      MobilityModel              const & mobility_;

    }; // current_density_on_edge

    template <typename DeviceType, typename SimulatorQuantity>
    class macroscopic_carrier_mask_filter
    {
      public:
        typedef bool    value_type;

        macroscopic_carrier_mask_filter(SimulatorQuantity const & quan) : quantity_(quan) {}

        value_type operator()(viennagrid_element_id c) const { return quantity_.get_unknown_mask(c); } //TODO: Might need fixing!

      private:
        SimulatorQuantity const & quantity_;
    };

    template <typename DeviceType, typename ValueT>
    class macroscopic_carrier_mask_filter<DeviceType, viennashe::const_quantity<ValueT> >
    {
      public:
        typedef bool    value_type;

        macroscopic_carrier_mask_filter(viennashe::const_quantity<ValueT> const & quan) : quantity_(quan) {}

        value_type operator()(viennagrid_element_id c) const { return true; } //TODO: Might need fixing!

      private:
        viennashe::const_quantity<ValueT> const & quantity_;
    };

  } // namespace detail


  /**
   * @brief An accessor to the current density on vertices and edges (drift diffusion only!)
   */
  template <typename DeviceType,
            typename PotentialQuantityType,
            typename CarrierQuantityType,
            typename MobilityModel>
  class current_density_wrapper
  {
    public:
      typedef std::vector<double> value_type;

      current_density_wrapper(DeviceType const & device,
                              viennashe::carrier_type_id ctype,
                              PotentialQuantityType const & potential,
                              CarrierQuantityType const & carrier,
                              MobilityModel const & mobility_model)
        : device_(device), carrier_type_id_(ctype), potential_(potential), carrier_(carrier), mobility_(mobility_model) { }

        value_type operator()(viennagrid_element_id cell_or_facet) const
        {
          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim));

          viennagrid_dimension element_dim = viennagrid_topological_dimension_from_element_id(cell_or_facet);

          if (cell_dim == element_dim) // cell provided
          {
            typedef detail::current_density_on_facet<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> current_density_evaluator;

            detail::macroscopic_carrier_mask_filter<DeviceType, CarrierQuantityType>    carrier_contribution_filter(carrier_);

            std::vector<double> J(3);

            if (!carrier_contribution_filter(cell_or_facet)) return J;

            viennashe::util::value_holder_functor<std::vector<double> > result;

            current_density_evaluator facet_evaluator(device_, carrier_type_id_, potential_, carrier_, mobility_);

            viennashe::util::dual_box_flux_to_cell(device_, cell_or_facet, result, facet_evaluator);

            return result();
          }
          else if (cell_dim == element_dim + 1) // facet provided
          {
            detail::current_density_on_facet<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> edge_evaluator(device_, carrier_type_id_, potential_, carrier_, mobility_);

            return edge_evaluator(cell_or_facet);
          }

          return value_type();
        }

    private:
      DeviceType                 const & device_;
      viennashe::carrier_type_id         carrier_type_id_;
      PotentialQuantityType      const & potential_;
      CarrierQuantityType        const & carrier_;
      MobilityModel              const & mobility_;
  };


  /** @brief Convenience function for writing the electric field to a container
   *
   * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param potential        An accessor to the electrostatic potential
   * @param carrier          The carrier quantity (n or p)
   * @param ctype            The carrier type identifier
   * @param mobility_model   The mobility model functor
   * @param container        A reference to the container, which is going to be filled with values, ie. std::vector
   */
  template <typename DeviceType,
            typename PotentialQuantityType,
            typename CarrierQuantityType,
            typename MobilityModel>
  void write_current_density_to_quantity_field(DeviceType const & device,
                                               PotentialQuantityType const & potential,
                                               CarrierQuantityType const & carrier,
                                               viennashe::carrier_type_id ctype,
                                               MobilityModel const & mobility_model,
                                               viennagrid_quantity_field field)
  {
    current_density_wrapper<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> Jfield(device, ctype, potential, carrier, mobility_model);

    viennashe::write_macroscopic_quantity_to_quantity_field(device, Jfield, field);
  }


  /**
   * @brief Checks current conservation for SHE. Writes information using log::info().
   * @param device             The device
   * @param current_on_facet   A function object returning the normal components of the current on each facet
   */
  template<typename DeviceT, typename CurrentDensityT>
  void check_current_conservation(DeviceT const & device,
                                  CurrentDensityT const & current_on_facet)
  {
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      double net_box_flux = 0;
      double partial_flux_max = 0;

      std::vector<viennagrid_numeric> centroid_cell(3);
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device.mesh(), *cit, &(centroid_cell[0])));

      viennagrid_element_id *facets_on_cell_begin, *facets_on_cell_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_boundary_elements(device.mesh(), *cit, cell_dim - 1, &facets_on_cell_begin, &facets_on_cell_end));

      for (viennagrid_element_id *focit  = facets_on_cell_begin;
                                  focit != facets_on_cell_end;
                                ++focit)
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


        double contribution_from_facet = current_on_facet(*focit)[0] * weighted_interface_area;
        if (std::fabs(contribution_from_facet) > partial_flux_max)
          partial_flux_max = std::fabs(contribution_from_facet);

        viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), *focit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

        if (*cit == cells_on_facet_begin[0]) //reference direction is into box
          net_box_flux -= contribution_from_facet;
        else
          net_box_flux += contribution_from_facet;
      }// for facets

      if (partial_flux_max > 0)
      {
        if ( std::fabs(net_box_flux / partial_flux_max) > 1e-10)
        {
          log::info() << "-- Box " << viennagrid_index_from_element_id(*cit) << " --" << std::endl;
          log::info() << "Net box flux (should be zero, unless at contact): " << net_box_flux << " with maximum " << partial_flux_max << std::endl;
        }
      }

    } // for cells

  }

  template <typename DeviceType,
            typename PotentialQuantityType,
            typename CarrierQuantityType,
            typename MobilityModel>
  void check_current_conservation(DeviceType const & device,
                                  viennashe::carrier_type_id ctype,
                                  PotentialQuantityType const & potential,
                                  CarrierQuantityType const & carrier,
                                  MobilityModel const & mobility_model)
  {
    current_density_wrapper<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> Jfield(device, ctype, potential, carrier, mobility_model);

    viennashe::check_current_conservation(device, Jfield);
  }


} //namespace viennashe

#endif
