#ifndef VIENNASHE_UTIL_POSTPROC_ELECTRIC_FIELD_HPP
#define VIENNASHE_UTIL_POSTPROC_ELECTRIC_FIELD_HPP

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
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/util/dual_box_flux.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"

/** @file viennashe/postproc/electric_field.hpp
    @brief Computes the electric field from a potential
 */

namespace viennashe
{
  namespace detail
  {
    /** @brief An accessor to the electric field along a given edge. Electrostatic potential required */
    /* TODO: Implement
    template <typename DeviceType,
              typename PotentialAccessorType>
    struct electric_field_on_facet
    {
      typedef typename DeviceType::mesh_type MeshType;

      typedef typename viennagrid::result_of::facet<MeshType>::type     FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

      electric_field_on_facet(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) {}

      template < typename FacetType>
          double operator()(FacetType const & facet) const
      {
        typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type    CellOnFacetContainer;
        typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

        CellOnFacetContainer cells_on_facet(device_.mesh(), viennagrid::handle(device_.mesh(), facet));

        if (cells_on_facet.size() < 2)
          return 0;

        CellOnFacetIterator cofit = cells_on_facet.begin();
        CellType const & c1 = *cofit;
        ++cofit;
        CellType const & c2 = *cofit;


        const double potential_center = potential_.get_value(c1);
        const double potential_outer  = potential_.get_value(c2);

        const double connection_len  = viennagrid::norm_2(viennagrid::centroid(c2) - viennagrid::centroid(c1));
        const double Emag     = -(potential_outer - potential_center) / connection_len;

        return Emag;
      }
      private:
        DeviceType const & device_;
        PotentialAccessorType const & potential_;

    }; // electric_field_on_edge
    */

  } // namespace detail

  /** @brief An accessor to the electric field on vertices and edges. Potential requiered */
  template < typename DeviceType, typename PotentialAccessorType >
  struct electric_field_wrapper
  {
  public:
    typedef std::vector<double> value_type;

    electric_field_wrapper(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) { }

    std::vector<double> operator()(viennagrid_element_id facet) const
    {
      //viennashe::detail::electric_field_on_facet<DeviceType, PotentialAccessorType> facet_eval(device_, potential_);

      throw std::runtime_error("electric_field_wrapper::operator(): TODO: Implement");

      std::vector<double> ret(3);
      return ret; //facet_eval(facet);
    }

/*
    value_type operator()(CellType const & cell) const
    {
      typedef typename viennagrid::result_of::const_facet_range<CellType>::type   FacetOnCellContainer;
      typedef detail::electric_field_on_facet<DeviceType, PotentialAccessorType>  FieldOnFacetEvaluator;

      viennashe::materials::checker no_conductor_filter(MATERIAL_NO_CONDUCTOR_ID);

      std::vector<double> E(3);

      if (!no_conductor_filter(device_.get_material(cell))) return E;

      viennashe::util::value_holder_functor<std::vector<double> > result;

      FieldOnFacetEvaluator facet_evaluator(device_, potential_);

      FacetOnCellContainer facets_on_cell(cell);
      viennashe::util::dual_box_flux_to_cell(device_,
                                             cell, facets_on_cell,
                                             result, facet_evaluator);

      return result();
    } // operator()
*/
    private:
      DeviceType const & device_;
      PotentialAccessorType const & potential_;


  }; // electric_field_wrapper


  /** @brief Convenience function for writing the electric field to a container
   *
   * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param potential        A quantity accessor to the potential
   * @param container        The container to be filled with values
   * @tparam ContainerType   A container type, for example: std::vector or std::deque
   */
  template <typename DeviceType,
            typename PotentialAccessor>
  void write_electric_field_to_quantity_field(DeviceType const & device,
                                              PotentialAccessor const & potential,
                                              viennagrid_quantity_field field)
  {
    electric_field_wrapper<DeviceType, PotentialAccessor> Efield(device, potential);

    viennashe::write_macroscopic_quantity_to_quantity_field(device, Efield, field);
  }

} // viennashe

#endif

