#ifndef VIENNASHE_UTIL_POSTPROC_DISPLACEMENT_CURRENT_HPP
#define VIENNASHE_UTIL_POSTPROC_DISPLACEMENT_CURRENT_HPP

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
#include "viennashe/accessors.hpp"

#include "viennashe/postproc/electric_field.hpp"

/** @file viennashe/postproc/displacement_current.hpp
    @brief Computes the electric field from a potential
 */

namespace viennashe
{
  namespace detail
  {
    // TODO: Revive

    /** @brief An accessor to the displacement current along edges */
    /*
    template < typename DeviceType, typename PotentialAccessorType, typename OldPotentialAccessorType >
    struct displacement_current_on_edge
    {

    private:
      typedef typename DeviceType::mesh_type  MeshType;

      typedef typename viennagrid::result_of::point<MeshType>::type   PointType;

    public:
      typedef typename viennagrid::result_of::vertex<MeshType>::type  VertexType;
      typedef typename viennagrid::result_of::edge<MeshType>::type    EdgeType;

      displacement_current_on_edge(DeviceType const & device,
                                   PotentialAccessorType const & potential,
                                   OldPotentialAccessorType const & old_potential,
                                   double rtime
                                  )
        : device_(device), potential_(potential), old_potential_(old_potential), rtime_(rtime)
      { }

      double operator()(EdgeType const & edge) const
      {
        typedef viennashe::util::no_contact_filter<DeviceType>  PotentialFilterType;

        if(rtime_ <= 0)
          return 0;

        PotentialFilterType                           potential_filter(device_);
        viennashe::permittivity_accessor<DeviceType>  permittivity(device_);

        viennashe::detail::electric_field_on_facet<DeviceType, PotentialAccessorType> efield_new(device_, potential_);
        viennashe::detail::electric_field_on_facet<DeviceType, PotentialAccessorType> efield_old(device_, old_potential_);

        const double A = util::effective_voronoi_interface(edge, potential_filter);

        if (A == 0) return 0; // safety check

        const double eps = util::effective_weighted_voronoi_interface(edge, potential_filter, permittivity) / A; // (sum_i A_i * eps_i) / (sum_i A_i) ... to get the average eps

        const double Dmag_new = efield_new(edge) * eps;
        const double Dmag_old = efield_old(edge) * eps;

        const double current_density = rtime_ * ( Dmag_new - Dmag_old );

        return current_density;
      } // operator()

    private:
      DeviceType               const & device_;
      PotentialAccessorType    const & potential_;
      OldPotentialAccessorType const & old_potential_;
      double                   rtime_;


    }; // displacement_current_on_edge
    */


  } // namespace detail

  /** @brief An accessor to the displacement current on edges and vertices */
  /*
  template < typename DeviceType, typename PotentialAccessorType, typename OldPotentialAccessorType >
  struct displacement_current_wrapper
  {
    private:
      typedef typename DeviceType::mesh_type  MeshType;

      typedef typename viennagrid::result_of::point<MeshType>::type    PointType;

    public:
      typedef typename viennagrid::result_of::vertex<MeshType> ::type  VertexType;
      typedef typename viennagrid::result_of::edge<MeshType> ::type    EdgeType;

      typedef std::vector<double> value_type;

      displacement_current_wrapper(DeviceType const & device,
                                   PotentialAccessorType const & potential,
                                   OldPotentialAccessorType const & old_potential,
                                   double rtime
                                  )
        : _device(device), _potential(potential), _old_potential(old_potential), _rtime(rtime)
      { }

      double operator()(EdgeType const & edge) const
      {
        viennashe::detail::displacement_current_on_edge<DeviceType, PotentialAccessorType, OldPotentialAccessorType> eval(_device, _potential, _old_potential, _rtime);
        return eval(edge);
      }

      value_type operator()(VertexType const & vertex) const
      {
        typedef typename viennagrid::result_of::const_edge_range<VertexType>::type      EdgeOnVertexContainer;

        typedef viennashe::util::no_contact_filter<DeviceType>  PotentialFilterType;

        PotentialFilterType  poisson_filter(_device);

        viennashe::detail::displacement_current_on_edge<DeviceType, PotentialAccessorType, OldPotentialAccessorType> edge_eval(_device, _potential, _old_potential, _rtime);

        std::vector<double> Jd(3);
        Jd[0] = 0; Jd[1] = 0; Jd[2] = 0;

        viennashe::util::value_holder_functor<PointType> result;
        EdgeOnVertexContainer edges_on_vertex(vertex, _device.mesh());
        viennashe::util::dual_box_flux_to_vertex(vertex, edges_on_vertex,
                                                 result, edge_eval,
                                                 poisson_filter, poisson_filter);

        for (std::size_t i=0; i < (PointType::dim) ; ++i) Jd[i] = result()[i];

        return Jd;

      } // operator()

    private:
      DeviceType               const & _device;
      PotentialAccessorType    const & _potential;
      OldPotentialAccessorType const & _old_potential;
      double                  _rtime;


  }; // displacement_current_wrapper
  */

  /** @brief Convenience function for writing the displacement current to the mesh.
   *
   * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
   * @param rtime            1.0/dt
   * @param potential        The quantity accessor for the potential of the current time step
   * @param old_potential    The quantity accessor for the potential of the old time step
   * @param quantity_name    String identifying the current density
   */
  /*
  template <typename DeviceType,
            typename PotentialAccessor,
            typename OldPotentialAccessorType,
            typename ContainerType>
  void write_displacement_current_to_container(DeviceType const & device,
                                               PotentialAccessor const & potential,
                                               PotentialAccessor const & old_potential,
                                               double rtime,
                                               ContainerType & container)
  {
    displacement_current_wrapper<DeviceType, PotentialAccessor, OldPotentialAccessorType> Jd(device, potential, old_potential, rtime);

    viennashe::she::write_macroscopic_quantity_to_container(device, Jd, container);
  }*/


} // viennashe

#endif

