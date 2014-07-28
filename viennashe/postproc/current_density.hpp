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
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/volume.hpp"
#include "viennagrid/algorithm/interface.hpp"

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
      typedef typename DeviceType::mesh_type MeshType;

      typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
      typedef typename viennagrid::result_of::vertex<MeshType>::type    VertexType;
      typedef typename viennagrid::result_of::facet<MeshType>::type     FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

      typedef double value_type;

      current_density_on_facet(DeviceType const & device,
                               viennashe::carrier_type_id ctype,
                               PotentialAccessorType const & potential,
                               AccessorTypeCarrier const & carrier,
                               MobilityModel const & mobility_model)
      : device_(device), carrier_type_id_(ctype), potential_(potential), carrier_(carrier), mobility_(mobility_model) { }

      value_type operator()(FacetType const & facet) const
      {
        typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type    CellOnFacetContainer;
        typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

        typename viennashe::contact_carrier_density_accessor<DeviceType> bnd_carrier_density(device_, carrier_type_id_);

        scharfetter_gummel flux_approximator(carrier_type_id_);

        CellOnFacetContainer cells_on_facet(device_.mesh(), viennagrid::handle(device_.mesh(), facet));

        if (cells_on_facet.size() < 2)
          return 0;

        CellOnFacetIterator cofit = cells_on_facet.begin();
        CellType const & c1 = *cofit;
        ++cofit;
        CellType const & c2 = *cofit;

        if (viennashe::materials::is_insulator(device_.get_material(c1)) || viennashe::materials::is_insulator(device_.get_material(c2))) // no current into insulator
          return 0;

        // at least one of the two cells must be a semiconductor:
        if (!viennashe::materials::is_semiconductor(device_.get_material(c1)) && !viennashe::materials::is_semiconductor(device_.get_material(c2)))
          return 0;

        const double potential_center = potential_(c1);
        const double carrier_center   = viennashe::materials::is_conductor(device_.get_material(c1))
                                        ? bnd_carrier_density(c2, device_.get_lattice_temperature(c2))
                                        : carrier_.get_value(c1);
        const double T  = device_.get_lattice_temperature(c1);

        const double mobility = mobility_.evaluate(c1, facet, c2);

        const double connection_len  = viennagrid::norm_2(viennagrid::centroid(c2) - viennagrid::centroid(c1));
        const double potential_outer = potential_(c2);
        const double carrier_outer   = viennashe::materials::is_conductor(device_.get_material(c2))
                                       ? bnd_carrier_density(c1, device_.get_lattice_temperature(c1))
                                       : carrier_.get_value(c2);

        if ( carrier_outer <= 0 || carrier_center <= 0 ) return 0;

        const double polarity = (carrier_type_id_ == viennashe::ELECTRON_TYPE_ID) ? -1.0 : 1.0;
        const double charge_flux = flux_approximator(carrier_center, carrier_outer, potential_center, potential_outer, connection_len, mobility, T);
        const double Jmag = polarity * mobility * charge_flux;

        return Jmag;
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
        typedef typename DeviceType::mesh_type           MeshType;

      public:
        typedef typename viennagrid::result_of::cell<MeshType>::type    cell_type;
        typedef bool    value_type;

        macroscopic_carrier_mask_filter(DeviceType const & dev, SimulatorQuantity const & quan) : device_(dev), quantity_(quan) {}

        value_type operator()(cell_type const & c) const
        {
          return viennashe::materials::is_semiconductor(device_.get_material(c));
          //return quantity_.get_unknown_mask(c); //TODO: Might need fixing!
        }

      private:
        DeviceType const & device_;
        SimulatorQuantity const & quantity_;
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
    protected:
      typedef typename DeviceType::mesh_type MeshType;

      typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
      typedef typename viennagrid::result_of::vertex<MeshType>::type    VertexType;
      typedef typename viennagrid::result_of::facet<MeshType>::type     FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

    public:
      typedef std::vector<double> value_type;

      current_density_wrapper(DeviceType const & device,
                              viennashe::carrier_type_id ctype,
                              PotentialQuantityType const & potential,
                              CarrierQuantityType const & carrier,
                              MobilityModel const & mobility_model)
        : device_(device), carrier_type_id_(ctype), potential_(potential), carrier_(carrier), mobility_(mobility_model) { }

        double operator()(FacetType const & facet) const
        {
          typedef detail::current_density_on_facet<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> current_density_evaluator;
          current_density_evaluator edge_evaluator(device_, carrier_type_id_, potential_, carrier_, mobility_);
          return edge_evaluator(facet);
        }

        value_type operator()(CellType const & cell) const
        {
          typedef typename viennagrid::result_of::const_facet_range<CellType>::type   FacetOnCellContainer;
          typedef detail::current_density_on_facet<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> current_density_evaluator;

          detail::macroscopic_carrier_mask_filter<DeviceType, CarrierQuantityType>    carrier_contribution_filter(device_, carrier_);

          std::vector<double> J(3);

          if (!carrier_contribution_filter(cell)) return J;

          viennashe::util::value_holder_functor<std::vector<double> > result;

          current_density_evaluator facet_evaluator(device_, carrier_type_id_, potential_, carrier_, mobility_);

          FacetOnCellContainer facets_on_cell(cell);
          viennashe::util::dual_box_flux_to_cell(device_,
                                                 cell, facets_on_cell,
                                                 result, facet_evaluator);

          return result();
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
  template <typename DeviceT,
            typename ContainerType>
  void write_current_density_to_container(viennashe::simulator<DeviceT> const & sim,
                                          viennashe::carrier_type_id ctype,
                                          ContainerType & container)
  {
    typedef typename viennashe::simulator<DeviceT>::potential_type         PotentialAccessor;
    typedef typename viennashe::simulator<DeviceT>::electron_density_type  CarrierAccessor;
    typedef viennashe::config::mobility_type                               MobilityModel;

    typedef viennashe::current_density_wrapper<DeviceT, PotentialAccessor, CarrierAccessor, MobilityModel>  CurrentDensityWrapper;
    CurrentDensityWrapper Jfield(sim.device(), ctype, sim.potential(), (ctype == ELECTRON_TYPE_ID) ? sim.electron_density() : sim.hole_density(), sim.config().mobility(ctype));

    viennashe::write_macroscopic_quantity_to_container(sim.device(), Jfield, container);
  }


  /**
   * @brief Checks current conservation for SHE. Writes information using log::info().
   * @param device             The device
   * @param current_on_facet   A function object returning the normal components of the current on each facet
   */
  template <typename DeviceT, typename CurrentDensityT>
  void check_current_conservation(DeviceT const & device,
                                  CurrentDensityT const & current_on_facet)
  {
    typedef typename DeviceT::mesh_type                 MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
    typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;
    typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

    typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;

    CellContainer cells(device.mesh());
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {
      double net_box_flux = 0;
      double partial_flux_max = 0;

      PointType centroid_cell = viennagrid::centroid(*cit);

      FacetOnCellContainer facets_on_cell(*cit);
      for (FacetOnCellIterator focit = facets_on_cell.begin();
          focit != facets_on_cell.end();
          ++focit)
      {
        CellType const *other_cell_ptr = util::get_other_cell_of_facet(device.mesh(), *focit, *cit);

        if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

        PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
        PointType cell_connection = centroid_other_cell - centroid_cell;
        PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
        PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);
        const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);
        double contribution_from_facet = current_on_facet(*focit) * weighted_interface_area;
        if (std::fabs(contribution_from_facet) > partial_flux_max)
          partial_flux_max = std::fabs(contribution_from_facet);

        CellOnFacetContainer cells_on_facet(device.mesh(), focit.handle());
        if ( &(*cit) == &(cells_on_facet[0]) ) //reference direction is into box
          net_box_flux -= contribution_from_facet;
        else
          net_box_flux += contribution_from_facet;
      }// for facets

      if (partial_flux_max > 0)
      {
        if ( std::fabs(net_box_flux / partial_flux_max) > 1e-10)
        {
          log::info() << "-- Box " << *cit << " --" << std::endl;
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
