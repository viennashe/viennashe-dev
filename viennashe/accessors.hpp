#ifndef VIENNASHE_ACCESSORS_HPP
#define VIENNASHE_ACCESSORS_HPP

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

// viennagrid
#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/mesh/neighbor_iteration.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/exception.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/device.hpp"

#include "viennashe/util/misc.hpp"

/** @file  accessors.hpp
    @brief Contains the definition of per-device accessors (read-only!) for various quantities.
*/


namespace viennashe
{
  namespace detail
  {
    /**
     * @brief Retrieves the averaged dopings from neighboring semiconductor cells
     * @param device The device
     * @param cell A cell in the device
     * @param doping_n_ret Return by reference; The averaged donor doping
     * @param doping_p_ret Return by reference; The averaged acceptor doping
     */
    template < typename DeviceType, typename CellType  >
    void get_dopings_from_neighboring_semiconductor_cells(DeviceType const & device, CellType const & cell,
                                                          double * doping_n_ret, double * doping_p_ret)
    {
      typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
      typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

      // Ensure unity doping by default
      double doping_n = 1;
      double doping_p = 1;

      long N = 0;
      FacetOnCellContainer facets(cell);
      for (FacetOnCellIterator focit = facets.begin(); focit != facets.end(); ++focit)
      {
        CellType const *other_cell_ptr = viennashe::util::get_other_cell_of_facet(device.mesh(), *focit, cell);
        if (other_cell_ptr == 0) continue;

        if ( viennashe::materials::is_semiconductor(device.get_material(*other_cell_ptr)) )
        {
          doping_n *= device.get_doping_n(*other_cell_ptr);
          doping_p *= device.get_doping_p(*other_cell_ptr);
          N += 1;
        }
      }
      if(N > 1) doping_n = std::pow(doping_n, 1.0/N);
      if(N > 1) doping_p = std::pow(doping_p, 1.0/N);

      // RETURN
      if(doping_n_ret) *doping_n_ret = doping_n;
      if(doping_p_ret) *doping_p_ret = doping_p;
    }

  } // namespace detail

  //
  // Accessors
  //
  /** @brief An accessor which returns a constant value independent of the object passed */
  template <typename ConstantType>
  class constant_accessor
  {
    public:
      typedef ConstantType    value_type;

      constant_accessor(ConstantType const & c) : constant_(c) {}

      template <typename T>
      value_type operator()(T const &) const { return constant_; }

    private:
      ConstantType constant_;
  };

  /** @brief Returns the lattice temperature on the device */
  template <typename DeviceType>
  class lattice_temperature_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;
    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      lattice_temperature_accessor(DeviceType const & d) : device_(d) {}

      value_type operator()(cell_type const & c) const { return device_.get_lattice_temperature(c); }

    private:
      DeviceType const & device_;
  };

  /** @brief Returns the thermal potential in the device */
  template <typename DeviceType>
  class thermal_potential_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;
    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      thermal_potential_accessor(DeviceType const & d) : device_(d) {}

      value_type operator()(cell_type const & c) const { return viennashe::physics::get_thermal_potential(device_.get_lattice_temperature(c)); }

    private:
      DeviceType const & device_;
  };

  /** @brief Accessor for returning the doping (donator/acceptor doping is defined in the CTOR) */
  template <typename DeviceType>
  class doping_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;
    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      /** @brief This CTOR specifies donator doping */
      doping_accessor(DeviceType const & d, viennashe::carrier_type_id ctype) : device_(d), is_doping_n_(ctype == viennashe::ELECTRON_TYPE_ID) {}

      value_type operator()(cell_type   const & c) const { return is_doping_n_ ? device_.get_doping_n(c) : device_.get_doping_p(c); }

    private:
      DeviceType const & device_;
      bool is_doping_n_;
  };


  /** @brief Accessor for retrieving the built-in potential in the device */
  template <typename DeviceType>
  class built_in_potential_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef typename viennagrid::result_of::facet<MeshType>::type     facet_type;
      typedef double    value_type;

      built_in_potential_accessor(DeviceType const & d) : device_(d) {}

      value_type operator()(cell_type const & c) const
      {
        if (viennashe::materials::is_semiconductor(device_.get_material(c)))
          return viennashe::physics::built_in_potential(device_.get_lattice_temperature(c),
                                                        device_.get_doping_n(c),
                                                        device_.get_doping_p(c));

        if (viennashe::materials::is_conductor(device_.get_material(c))) // important for boundary conditions
        {
          typedef typename viennagrid::result_of::const_neighbor_range<MeshType, cell_type, facet_type>::type    NeighborRange;
          typedef typename viennagrid::result_of::iterator<NeighborRange>::type                                  NeighborIterator;

          // for a cell attached to a semiconductor, return the contact potential corrected by the built-in potential:
          NeighborRange neighbors(device_.mesh(), viennagrid::handle(device_.mesh(), c));
          for (NeighborIterator nit = neighbors.begin(); nit != neighbors.end(); ++nit)
          {
            if (viennashe::materials::is_semiconductor(device_.get_material(*nit)))
              return viennashe::physics::built_in_potential(device_.get_lattice_temperature(*nit),
                                                            device_.get_doping_n(*nit),
                                                            device_.get_doping_p(*nit));
          }
        }

        // cell not attached to semiconductor, so we can safely assume no built-in potential
        return 0.0;
      }

    private:
      DeviceType const & device_;
  };


  /** @brief Accessor for obtaining the Dirichlet boundary condition for the Poisson equation in the device (contact-potential plus built-in potential) */
  template <typename DeviceT>
  class boundary_potential_accessor
  {
      typedef typename DeviceT::mesh_type     MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      boundary_potential_accessor(DeviceT const & d) : device_(d), built_in_pot_(d) {}

      value_type operator()(cell_type const & c) const
      {
        if (viennashe::materials::is_semiconductor(device_.get_material(c))
            ||  viennashe::materials::is_conductor(device_.get_material(c)))
          return device_.get_contact_potential(c) + built_in_pot_(c);

        throw std::runtime_error("Logic error: Accessing boundary potential for insulators!");
      }

    private:
      DeviceT const & device_;
      built_in_potential_accessor<DeviceT> built_in_pot_;
  };

  /** @brief Accessor for obtaining the contact potential in the device */
  template <typename DeviceType>
  class contact_potential_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      contact_potential_accessor(DeviceType const & d) : device_(d) {}

      value_type operator()(cell_type const & c) const
      {
        return device_.get_contact_potential(c);
      }

    private:
      DeviceType const & device_;
  };


  /** @brief Accessor for obtaining the permittivity in the device */
  template <typename DeviceType>
  class permittivity_accessor
  {
      typedef typename DeviceType::mesh_type     MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double    value_type;

      permittivity_accessor(DeviceType const & d) : device_(d) {}

      value_type operator()(cell_type const & cell) const { return viennashe::materials::permittivity(device_.get_material(cell)); }

    private:
      DeviceType const & device_;
  };


  /** @brief Accessor for fixed charges */
  template <typename DeviceType>
  class fixed_charge_accessor
  {
      typedef typename DeviceType::mesh_type    MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
      typedef double  value_type;

      fixed_charge_accessor(DeviceType const & d) : device_(d) { }

      /** @brief Returns the fixed charge in As */
      value_type operator()(cell_type const & c) const
      {
        return device_.get_fixed_charge(c);
      }

    private:
     DeviceType const & device_;
  };

  /** @brief Accessor to get the diffusivity. Used in the assembly of the heat diffusion equation */
  template <typename DeviceType>
  class diffusivity_accessor
  {
    typedef typename DeviceType::mesh_type     MeshType;

  public:
    typedef typename viennagrid::result_of::cell<MeshType>::type      cell_type;
    typedef double    value_type;

    diffusivity_accessor(DeviceType const & d) : device_(d) {}

    value_type operator()(cell_type const & cell, double T) const
    {
      (void)T; // unused at the moment
      return viennashe::materials::diffusivity(device_.get_material(cell));
    }

  private:
    DeviceType const & device_;
 };

 /** @brief Returns the carrier density at contacts modelled as thermal baths (used by DD and SHE) */
  template <typename DeviceType>
  struct contact_carrier_density_accessor
  {
    typedef double value_type;

    contact_carrier_density_accessor(DeviceType const & d, viennashe::carrier_type_id ctype, bool doping_bnd = false)
     : device_(d), carrier_type_id_(ctype), doping_bnd_(doping_bnd) { }

    template < typename ElementType >
    value_type operator()(ElementType const & el) const
    {
      const double temperature = device_.get_lattice_temperature(el);
      return (this->operator ()(el, temperature));
    }

    template < typename ElementType >
    value_type operator()(ElementType const & el, double T) const
    {
      if (T <= 0.0) throw viennashe::invalid_value_exception("contact_carrier_density_accessor: T <= 0!");

      const long matid = device_.get_material(el);

      if ( viennashe::materials::is_insulator(matid))
        return 0;

      double doping_n    = device_.get_doping_n(el);
      double doping_p    = device_.get_doping_p(el);
      if( viennashe::materials::is_conductor(matid) )
      {
        // ATTENTION: return by reference!
        viennashe::detail::get_dopings_from_neighboring_semiconductor_cells(device_, el, &doping_n, &doping_p);
      }

      if (doping_bnd_ == true)
      {
        if(carrier_type_id_ == viennashe::ELECTRON_TYPE_ID)
          return doping_n;
        else if(carrier_type_id_ == viennashe::HOLE_TYPE_ID)
          return doping_p;
        else
          throw viennashe::carrier_type_not_supported_exception("contact_carrier_density_accessor: ctype");
      }
      else
      {
        return viennashe::physics::contact_carrier_ohm(T, doping_n, doping_p, carrier_type_id_);
      }
    }

    private:
      DeviceType const & device_;
      viennashe::carrier_type_id carrier_type_id_;
      bool doping_bnd_;
  };


} // namespace viennashe

#endif

