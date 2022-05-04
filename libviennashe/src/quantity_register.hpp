#ifndef LIBVIENNASHE_QUANTITY_REGISTER_HPP
#define	LIBVIENNASHE_QUANTITY_REGISTER_HPP

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



/** @file libviennashe/src/quantity_register.hpp
    @brief Common routines for the registry of simulator dependent quantities
*/


// C++ includes
#include "libviennashe/src/viennashe_all.hpp"
#include "libviennashe/src/quantity_wrappers.hpp"

// C includes
#include "libviennashe/include/quantity.h"

#include <cstring>

namespace libviennashe
{
  namespace detail
  {

    /**
     * @brief Registers a current density from a DD simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param potential An accessor to the potential
     * @param carrier An accessor to the carrier density (n or p)
     * @param ctype The carrier type. Must match the carrier density carrier type
     * @param mobility_model A suitable mobility model (DD mobility)
     * @param name A unique name for the new quantity
     */
    template <typename DeviceType,
      typename PotentialQuantityType,
      typename CarrierQuantityType,
      typename MobilityModel>
    void register_DD_current_density(quan_register_internal & reg,
                                     DeviceType const & device,
                                     PotentialQuantityType const & potential,
                                     CarrierQuantityType const & carrier,
                                     viennashe::carrier_type_id ctype,
                                     MobilityModel const & mobility_model, std::string name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennashe::current_density_wrapper<DeviceType, PotentialQuantityType, CarrierQuantityType, MobilityModel> CurrentDensityAccessorType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;

      CurrentDensityAccessorType Jfield(device, ctype, potential, carrier, mobility_model);

      libviennashe::quantity::accessor_based_array_quantity_wrapper<DeviceType, CurrentDensityAccessorType, CellTag> quan_cell(Jfield, device, name);
      reg.cell_based.register_quan(quan_cell);
    }

    /**
     * @brief Registers a current density from a SHE simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param shequan An accessor to the SHE quantities
     * @param conf The simulator configuration
     * @param name A unique name for the new quantity
     */
    template <typename DeviceType, typename SHEQuanT, typename ConfigT>
      void register_SHE_current_density(quan_register_internal & reg,
                                        DeviceType const & device,
                                        SHEQuanT const & shequan,
                                        ConfigT const & conf, std::string name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennashe::she::current_density_wrapper<DeviceType, SHEQuanT> CurrentDensityAccessorType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;

      CurrentDensityAccessorType Jfield(device, conf, shequan);

      libviennashe::quantity::accessor_based_array_quantity_wrapper<DeviceType, CurrentDensityAccessorType, CellTag> quan_cell(Jfield, device, name);
      reg.cell_based.register_quan(quan_cell);
    }

    /**
     * @brief Registers the average carrier velocity for a SHE simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param shequan An accessor to the SHE quantities
     * @param conf The configuration of the simulator used for computing the quantity shequan
     * @param name A unique name for the new quantity
     */
    template <typename DeviceType, typename SHEQuanT, typename ConfigT>
    void register_average_velocity(quan_register_internal & reg,
                                   DeviceType const & device,
                                   SHEQuanT const & shequan,
                                   ConfigT const & conf, std::string name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennashe::she::carrier_velocity_wrapper<DeviceType, SHEQuanT> AccessorType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;

      AccessorType acc(device, conf, shequan);
      libviennashe::quantity::accessor_based_array_quantity_wrapper<DeviceType, AccessorType, CellTag> quan_cell(acc, device, name);
      reg.cell_based.register_quan(quan_cell);
    }

    /**
     * @brief Registers the average kinetic energy for a SHE simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param shequan An accessor to the SHE quantities
     * @param conf The simulator configuration
     * @param name A unique name for the new quantity
     */
    template <typename DeviceType, typename SHEQuanT, typename ConfigT>
    void register_average_energy(quan_register_internal & reg,
                                 DeviceType const & device,
                                 SHEQuanT const & shequan,
                                 ConfigT const & conf, std::string name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;
      typedef typename viennashe::she::carrier_energy_wrapper<SHEQuanT> AccessorType;

      AccessorType acc(conf, shequan);

      libviennashe::quantity::accessor_based_quantity_wrapper<DeviceType, AccessorType, CellTag> quan_cell(acc, device, name);
      reg.cell_based.register_quan(quan_cell);
    }

    /**
     * @brief Registers the average expansion order (SHE simulator)
     * @param reg The internal quantity register
     * @param device The device
     * @param shequan An accessor to the SHE quantities
     * @param name A unique name for the new quantity
     */
    template <typename DeviceType, typename SHEQuanT>
    void register_average_expansion_order(quan_register_internal & reg,
                                 DeviceType const & device,
                                 SHEQuanT const & shequan,
                                 std::string name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;
      typedef typename viennashe::she::average_expansion_order_wrapper<DeviceType, SHEQuanT> AccessorType;

      AccessorType acc(shequan, 0.0, 1.0);

      libviennashe::quantity::accessor_based_quantity_wrapper<DeviceType, AccessorType, CellTag> quan_cell(acc, device, name);
      reg.cell_based.register_quan(quan_cell);
    }



    /**
     * @brief Registers the electric field (vector) for any simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param potential An accessor to the electrostatic potential
     * @param quantity_name A unique name for the new quantity
     */
    template <typename DeviceType, typename PotentialAccessor>
    void register_electric_field(quan_register_internal & reg, DeviceType const & device,
                                 PotentialAccessor const & potential, std::string quantity_name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;
      typedef typename viennashe::electric_field_wrapper<DeviceType, PotentialAccessor> ElectricFieldAccessorType;
      ElectricFieldAccessorType Efield(device, potential);

      libviennashe::quantity::accessor_based_array_quantity_wrapper<DeviceType, ElectricFieldAccessorType, CellTag> quan_vertex(Efield, device, quantity_name);
      reg.cell_based.register_quan(quan_vertex);
    }

    /**
     * @brief Registers the electric flux density (vector) for any simulator
     * @param reg The internal quantity register
     * @param device The device
     * @param potential An accessor to the elec. potential
     * @param quantity_name A unique name for the new quantity
     */
    template <typename DeviceType, typename PotentialAccessor>
    void register_electric_flux(quan_register_internal & reg, DeviceType const & device,
      PotentialAccessor const & potential, std::string quantity_name)
    {
      typedef typename DeviceType::mesh_type  MeshType;
      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;
      typedef typename viennashe::electric_flux_wrapper<DeviceType, PotentialAccessor> AccessorType;
      typedef typename libviennashe::quantity::accessor_based_array_quantity_wrapper<DeviceType, AccessorType, CellTag> WrapperType;

      AccessorType Dfield(device, potential);
      WrapperType quan_cell(Dfield, device, quantity_name);
      reg.cell_based.register_quan(quan_cell);
    }

  } // namespace detail

  /**
   * @brief Main quantity regsiter function. Add your simulator quantities here
   * @param sim The simulator from which quantities can be obtained
   * @param reg The internal quantity registry
   */
  template < typename SimulatorT >
  void register_quans(SimulatorT const & sim, quan_register_internal & reg)
  {
    typedef typename SimulatorT::device_type         DeviceType;
    typedef typename SimulatorT::ResultQuantityType  ResultQuantityType;

    typedef typename DeviceType::mesh_type  MeshType;
    typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;

    DeviceType const & device = sim.device();

    //
    // Register device quantities
    //

    // Doping
    viennashe::doping_accessor<DeviceType> doping_n(device, viennashe::ELECTRON_TYPE_ID);
    viennashe::doping_accessor<DeviceType> doping_p(device, viennashe::HOLE_TYPE_ID);
    quantity::accessor_based_quantity_wrapper<DeviceType, viennashe::doping_accessor<DeviceType>, CellTag> doping_n_vertex(doping_n, device, "Donor doping concentration");
    quantity::accessor_based_quantity_wrapper<DeviceType, viennashe::doping_accessor<DeviceType>, CellTag> doping_p_vertex(doping_p, device, "Acceptor doping concentration");
    reg.cell_based.register_quan(doping_n_vertex);
    reg.cell_based.register_quan(doping_p_vertex);

    // Built-in Potential
    viennashe::built_in_potential_accessor<DeviceType> builtinpot(device);
    quantity::accessor_based_quantity_wrapper<DeviceType, viennashe::built_in_potential_accessor<DeviceType>, CellTag> builtinpot_vt(builtinpot, device, "Built-In potential");
    reg.cell_based.register_quan(builtinpot_vt);

    //
    // Register basic quantities
    //

    typedef typename quantity::accessor_based_quantity_wrapper<DeviceType, ResultQuantityType, CellTag> ResultQuantityWrapperType;

    // Potential
    ResultQuantityWrapperType pot_vertex(sim.potential(), device, viennashe::quantity::potential());
    reg.cell_based.register_quan(pot_vertex);

    // Electron Concentration
    ResultQuantityWrapperType n_vertex(sim.electron_density(), device, viennashe::quantity::electron_density());
    reg.cell_based.register_quan(n_vertex);

    // Hole Concentration
    ResultQuantityWrapperType p_vertex(sim.hole_density(), device, viennashe::quantity::hole_density());
    reg.cell_based.register_quan(p_vertex);

    // DG Electron Correction Potential
    ResultQuantityWrapperType dgn_vertex(sim.dg_pot_n(), device, viennashe::quantity::density_gradient_electron_correction());
    reg.cell_based.register_quan(dgn_vertex);

    // DG Hole Correction Potential
    ResultQuantityWrapperType dgp_vertex(sim.dg_pot_p(), device, viennashe::quantity::density_gradient_hole_correction());
    reg.cell_based.register_quan(dgp_vertex);

    // Lattice Temperature
    ResultQuantityWrapperType hde_lat_temp(sim.quantities().lattice_temperature(), device, viennashe::quantity::lattice_temperature());
    reg.cell_based.register_quan(hde_lat_temp);

    viennashe::lattice_temperature_accessor<DeviceType> lattice_temp(device);
    quantity::accessor_based_quantity_wrapper<DeviceType, viennashe::lattice_temperature_accessor<DeviceType>, CellTag>
        tl_dev(lattice_temp, device, "Transport lattice temperature");
    reg.cell_based.register_quan(tl_dev);

    //
    // Register postprocessing quantities
    //

    // Electric Field (V/m)
    detail::register_electric_field(reg, device, sim.potential(), "Electric field");
    // Electric Flux (As/m^2)
    detail::register_electric_flux(reg, device, sim.potential(), "Electric flux density");

    if (sim.config().with_electrons())
    {
      if (sim.config().get_electron_equation() == viennashe::EQUATION_CONTINUITY)
      {
        // Current Density Electrons
        detail::register_DD_current_density(reg, device,
          sim.quantities().potential(),
          sim.quantities().get_unknown_quantity(viennashe::quantity::electron_density()),
          viennashe::ELECTRON_TYPE_ID,
          viennashe::models::create_constant_mobility_model(device, 0.1430), "Electron current density DD");
      }
      else if (sim.config().get_electron_equation() == viennashe::EQUATION_SHE)
      {
        // Current Density Electrons
        detail::register_SHE_current_density(reg, device, sim.quantities().electron_distribution_function(),
                                             sim.config(), "Electron current density SHE");
        // Average Carrier Energy
        detail::register_average_energy(reg, device, sim.quantities().electron_distribution_function(),
                                        sim.config(), "Average electron energy SHE");
        // Average Carrier Velocity
        detail::register_average_velocity(reg, device, sim.quantities().electron_distribution_function(),
                                          sim.config(), "Average electron velocity SHE");

        // Average Expansion Order
        detail::register_average_expansion_order(reg, device, sim.quantities().electron_distribution_function(), "Average electron expansion order SHE");

      }
      else
      {
        // TODO: throw ?
        viennashe::log::warn() << "libviennashe: register_quans(): Unknown equation type!" << std::endl;
      }
    }

    if (sim.config().with_holes())
    {
      if (sim.config().get_hole_equation() == viennashe::EQUATION_CONTINUITY)
      {
        // Current Density Holes
        detail::register_DD_current_density(reg, device,
          sim.quantities().potential(),
          sim.quantities().get_unknown_quantity(viennashe::quantity::hole_density()),
          viennashe::HOLE_TYPE_ID,
          viennashe::models::create_constant_mobility_model(device, 0.0460), "Hole current density DD");
      }
      else if (sim.config().get_hole_equation() == viennashe::EQUATION_SHE)
      {
        // Current Density Holes
        detail::register_SHE_current_density(reg, device, sim.quantities().hole_distribution_function(),
                                             sim.config(), "Hole current density SHE");
        // Average Carrier Energy
        detail::register_average_energy(reg, device, sim.quantities().hole_distribution_function(),
                                        sim.config(), "Average hole energy SHE");
        // Average Carrier Velocity
        detail::register_average_velocity(reg, device, sim.quantities().hole_distribution_function(),
                                          sim.config(), "Average hole velocity SHE");
        // Average Expansion Order
        detail::register_average_expansion_order(reg, device, sim.quantities().hole_distribution_function(), "Average hole expansion order SHE");

      }
      else
      {
        // TODO: throw ?
        viennashe::log::warn() << "libviennashe: register_quans(): Unknown equation type!" << std::endl;
      }
    }


  } // register_quans


  /**
   * @brief Fills the given C-arrays with the EDF at a vertex
   * @param quan The SHE quantities
   * @param edfacc A energy distribution function wrapper for quan
   * @param cell The cell for which to extract the EDF
   * @param ekin Return value: single array. Will hold the kinetic energies
   * @param edf Return value: single array. Will hold the EDF values
   * @param len Return value: single value. Will hold the length of each array
   */
  template <typename SHEQuanT, typename DFWrapperT, typename CellType >
  void she_fill_edf_at_cell(SHEQuanT const & quan, DFWrapperT const & edfacc, CellType const & cell, double ** ekin, double ** edf, viennashe_index_type * len)
  {
    const std::size_t idx = std::size_t(cell.id().get());
    len[idx] = static_cast<viennashe_index_type>(quan.get_value_H_size() - 1); // set length
    for (std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H)
    {
      const double energy = quan.get_kinetic_energy(cell, index_H);
      const double f = edfacc(cell, energy, index_H);
      ekin[idx][index_H - 1] = energy;
      edf[idx][index_H - 1]  = f;
    }
  }


  /**
   * @brief Fills the given C-arrays ekin and dos with the DOS at a vertex
   * @param num The number of energies for which to fill
   * @param deltaeps The delta in kinetic energy
   * @param dosacc A C++ accessor to the DOS
   * @param vt The vertex at which to extract the DOS
   * @param ekin Return value: single array. Will hold the energy values
   * @param dos Return value: single array. Will hold the DOS values
   * @param len Return value. Will hold the length of ekin and dos
   */
  template <typename DOSAccessorT, typename VertexType >
  void she_fill_dos_at_cell(viennashe_index_type num, double deltaeps, DOSAccessorT const & dosacc, VertexType const & vt, double ** ekin, double ** dos, viennashe_index_type * len)
  {
    const std::size_t idx = std::size_t(vt.id().get());
    len[idx] = num;
    for (std::size_t index_H = 0; index_H < num; ++index_H)
    {
      const double energy = deltaeps * index_H;
      ekin[idx][index_H - 1] = energy;
      dos[idx][index_H - 1] = dosacc.density_of_states(energy);
    }
  }

  /**
   * @brief Fills the given C-arrays with the complete EDF
   * @param sim The SHE simulator (run will NOT be called!)
   * @param ctype The carrier type for which the EDF shall be returned
   * @param ekin Return value: single array. Will hold the kinetic energies
   * @param edf Return value: single array. Will hold the EDF values
   * @param len Return value: single value. Will hold the length of each array
   */
  template <typename SimulatorT >
  void she_fill_edf(SimulatorT const & sim, viennashe::carrier_type_id ctype, double ** ekin, double ** edf, viennashe_index_type * len)
  {
    typedef typename SimulatorT::device_type  DeviceType;
    typedef typename DeviceType::mesh_type    MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

    DeviceType const & device = sim.device();
    CellContainer cells(device.mesh());

    if (ctype == viennashe::ELECTRON_TYPE_ID)
    {
      for (CellIterator cit = cells.begin();
           cit != cells.end();
           ++cit)
      {
        she_fill_edf_at_cell(sim.quantities().electron_distribution_function(), sim.edf(ctype), *cit, ekin, edf, len);
      }
    }
    else
    {
      for (CellIterator cit = cells.begin();
           cit != cells.end();
           ++cit)
      {
        she_fill_edf_at_cell(sim.quantities().hole_distribution_function(), sim.edf(ctype), *cit, ekin, edf, len);
      }
    }
  } // she_fill_edf

  /**
   * @brief Fills the given C-arrays ekin and dos with the DOS for every vertex
   * @param sim The SHE simulator (run will NOT be called!)
   * @param ctype The carrier type for which the EDF shall be returned
   * @param ekin Return value: single array. Will hold the energy values
   * @param dos Return value: single array. Will hold the DOS values
   * @param len Return value. Will hold the length of ekin and dos
   */
  template < typename SimulatorT >
  void she_fill_dos(SimulatorT const & sim, viennashe::carrier_type_id ctype, double ** ekin, double ** dos, viennashe_index_type * len)
  {
    typedef typename SimulatorT::device_type  DeviceType;
    typedef typename DeviceType::mesh_type    MeshType;

    typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type     VertexContainer;
    typedef typename viennagrid::result_of::iterator<VertexContainer>::type        VertexIterator;

    DeviceType const & device = sim.device();
    VertexContainer vertices(device.mesh());

    const double ekinmax = viennashe::physics::convert::eV_to_joule(10.0); // TODO: Think about better maximum for the kinetic energy ...
    const viennashe_index_type  num = static_cast<viennashe_index_type>(std::ceil(ekinmax / sim.config().energy_spacing()));

    for (VertexIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
    {
      she_fill_dos_at_cell(num, sim.config().energy_spacing(), sim.config().dispersion_relation(ctype), *vit, ekin, dos, len);
    }

  } // she_fill_dos

  /**
   * @brief Fills the given C-arrays ekin and vg with the group velocity for every vertex
   * @param sim The SHE simulator (run will NOT be called!)
   * @param ctype The carrier type for which the EDF shall be returned
   * @param ekin Return value: single array. Will hold the energy values
   * @param vg Return value: single array. Will hold the group velocity values
   * @param len Return value. Will hold the length of ekin and dos
   */
  template < typename SimulatorT >
  void she_fill_group_velocity(SimulatorT const & sim, viennashe::carrier_type_id ctype, double ** ekin, double ** vg, viennashe_index_type * len)
  {
    typedef typename SimulatorT::device_type  DeviceType;
    typedef typename DeviceType::mesh_type    MeshType;

    typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type     VertexContainer;
    typedef typename viennagrid::result_of::iterator<VertexContainer>::type        VertexIterator;

    DeviceType const & device = sim.device();
    VertexContainer vertices(device.mesh());

    const double ekinmax  = viennashe::physics::convert::eV_to_joule(10.0); // TODO: Think about better maximum for the kinetic energy ...

    const viennashe_index_type num = static_cast<viennashe_index_type>(std::ceil(ekinmax / sim.config().energy_spacing()));

    for (VertexIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
    {
      const std::size_t idx = std::size_t(vit->id().get());
      len[idx] = num;
      for (std::size_t index_H = 0; index_H < num; ++index_H)
      {
        const double energy = sim.config().energy_spacing() * index_H;
        ekin[idx][index_H - 1] = energy;
        vg[idx][index_H - 1] =  sim.config().dispersion_relation(ctype).velocity(energy);
      }
    }

  } // she_fill_group_velocity


} // namespace libviennashe


#endif	/* LIBVIENNASHE_QUANTITY_REGISTER_HPP */

