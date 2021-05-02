#ifndef VIENNASHE_UTIL_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP
#define VIENNASHE_UTIL_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP

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
#include "viennagrid/algorithm/inner_prod.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/util/checks.hpp"

#include "viennashe/log/log.hpp"

#include "viennashe/util/misc.hpp"
#include "viennashe/util/dual_box_flux.hpp"

/** @brief Computes the electric flux density from a potential
 * @file viennashe/postproc/electric_flux_density.hpp
 */

namespace viennashe
{
  namespace detail
  {
    /**
     * @brief Simple accessor to get the electric flux density along an edge
     */
    template < typename DeviceType, typename PotentialAccessorType >
    struct electric_flux_on_facet
    {
      typedef typename DeviceType::mesh_type MeshType;

      typedef typename viennagrid::result_of::facet<MeshType>::type     FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

      electric_flux_on_facet(DeviceType const & device, PotentialAccessorType const & potential) : device_(device), potential_(potential) {}

      template <typename FacetType>
      double operator()(FacetType const & facet) const
      {
        typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type    CellOnFacetContainer;
        typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;
        typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;

        CellOnFacetContainer cells_on_facet(device_.mesh(), viennagrid::handle(device_.mesh(), facet));

        if (cells_on_facet.size() < 2)
          return 0;

        CellOnFacetIterator cofit = cells_on_facet.begin();
        CellType const & c1 = *cofit;
        ++cofit;
        CellType const & c2 = *cofit;


        viennashe::permittivity_accessor<DeviceType>  permittivity(device_);

        const double potential_center = potential_.get_value(c1);
        const double potential_outer  = potential_.get_value(c2);

        PointType centroid_other_cell = viennagrid::centroid(c2);
        PointType centroid_cell = viennagrid::centroid(c1);
        PointType cell_connection = centroid_other_cell - centroid_cell;
        //PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
        //PointType facet_normal = cell_connection_normalized; // TODO: Substitute facet normal calculation here

        const double connection_len = viennagrid::norm_2(cell_connection);
        //const double interface_area = viennagrid::volume(facet) * viennagrid::inner_prod(facet_normal, cell_connection_normalized);

        PointType facet_center = viennagrid::centroid(facet);  //TODO: Use intersection of facet plane with connection
        double connection_in_cell = viennagrid::norm(facet_center - centroid_cell);
        double connection_in_other_cell = viennagrid::norm(facet_center - centroid_other_cell);
        const double permittivity_mean = (connection_in_cell + connection_in_other_cell) /
                                           (connection_in_cell/permittivity(c1) + connection_in_other_cell/permittivity(c2));

        const double Emag  = -(potential_outer - potential_center) / connection_len;

        return Emag * permittivity_mean;
      }

      private:
        DeviceType            const & device_;
        PotentialAccessorType const & potential_;


    }; // electric_flux_on_facet

  } // namespace detail


  /**
   * @brief An electric flux accessor. Provides the electric flux along edges and on vertices given the electrostatic potential.
   */
  template < typename DeviceType, typename PotentialAccessorType >
  struct electric_flux_wrapper
  {
  private:
    typedef typename DeviceType::mesh_type   MeshType;
    typedef typename MeshType::config_type   ConfigType;

    typedef typename viennagrid::result_of::point<MeshType>::type PointType;

  public:
    typedef typename viennagrid::result_of::facet<MeshType>::type     FacetType;
    typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

    typedef std::vector<double> value_type;

    electric_flux_wrapper(DeviceType const & device, PotentialAccessorType const & potential)
      : device_(device), potential_(potential)
    { }

    double operator()(FacetType const & facet) const
    {
      viennashe::detail::electric_flux_on_facet<DeviceType, PotentialAccessorType> facet_eval(device_, potential_);

      return facet_eval(facet);
    }


    value_type operator()(CellType const & cell) const
    {
      typedef typename viennagrid::result_of::const_facet_range<CellType>::type   FacetOnCellContainer;
      typedef detail::electric_flux_on_facet<DeviceType, PotentialAccessorType>  FieldOnFacetEvaluator;

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

    private:
      DeviceType            const & device_;
      PotentialAccessorType const & potential_;

  }; // electric_flux_wrapper


} // viennashe

#endif /* VIENNASHE_DD_POSTPROC_ELECTRIC_FLUX_DENSITY_HPP */

