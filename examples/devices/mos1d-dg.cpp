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

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include <iostream>
#include <cstdlib>
#include <vector>

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid mesh configurations:
#include "viennagrid/viennagrid.h"


/** \example mos1d-dg.cpp

  This example covers the simulation of a purely 1d MOS (metal-SiO2-Si-metal) using the density gradient model for quantum corrections.
  The device structure is follows: \code

    -------------------------------------
    |       |       |  |        |       |
    | Metal | Oxide |Si|   Si   | Metal |
    |       |       |  |        |       |
    -------------------------------------
    #   0       1    2     3       4        (segment IDs)

  \endcode

  Note that the silicon-layer is split into a 20nm channel with smaller grid spacing and a bulk part with larger grid spacing.

<h2>First Step: Generate a Mesh for the Device</h2>

To generate a suitable mesh for our MOS device, we use the built-in 1d/2d mesh generator in ViennaSHE.
In this example we populate the mesh generator configuration accordingly and then pass it to the device,
which then generates the mesh automatically.
**/

template <typename DeviceT>
void generate_mos1d_mesh(DeviceT & mos_device,
                         double len_gate,  // gate length
                         double cs_gate,   // cell size in the gate
                         double len_oxide, // oxide length
                         double cs_ox,     // cell size in the oxide
                         double len_bulk,  // length of the silicon bulk
                         double cs_bulk)   // cell size in the silicon bulk
{
  viennashe::util::device_generation_config gconf;

  // The gate contact is segment #0:
  gconf.add_segment(0,                                                       // starting location (meter)
                    len_gate,                                                // length (meter)
                    static_cast<unsigned long>(len_gate/cs_gate + 1.0) );    // number of mesh points in the segment

  // The oxide is segment #1
  gconf.add_segment(len_gate,
                    len_oxide,
                    static_cast<unsigned long>(len_oxide/cs_ox)  );

  // The channel is segment #2 and has a fixed length of 20 nm and a cell size of 0.05 nm
  gconf.add_segment(len_gate + len_oxide,
                    20e-9, // 20nm
                    static_cast<unsigned long>(+20e-9/0.05e-9) );

  // The bulk silicon is segment #3
  gconf.add_segment(len_gate + len_oxide + 20e-9,
                    len_bulk,
                    static_cast<unsigned long>(len_bulk/cs_bulk) );

  // The last segment #4 is the bulk contact
  gconf.add_segment(len_gate + len_oxide + 20e-9 + len_bulk,
                    2e-9,
                    static_cast<unsigned long>(2e-9/cs_bulk + 1.0) );

  // Pass the
  mos_device.generate_mesh(gconf);
}


/**
<h2>Second Step: Initialize the Device</h2>

Once the device mesh is generated, we need to initialize the various
device segments (aka. submeshes) with material parameters, contact
voltages, etc. For simplicity, we collect this initialization
in a separate function:
**/
template <typename DeviceType>
void init_device(DeviceType & device,
                 double Vg_init,  // The initial gate voltage
                 double Nd,       // The donor doping (m^-3)
                 double Na)       // The acceptor doping (m^-3)
{
  typedef typename DeviceType::segment_type        SegmentType;

  // Get the segments of the mesh
  SegmentType const & gate     = device.segment(0);
  SegmentType const & oxide    = device.segment(1);
  SegmentType const & silicon  = device.segment(2);
  SegmentType const & silicon2 = device.segment(3);
  SegmentType const & bulk     = device.segment(4);

  std::cout << "* init_device(): Setting material ..." << std::endl;

  /**
    First we set the materials per segment
  **/
  device.set_material(viennashe::materials::sio2(),  oxide);
  device.set_material(viennashe::materials::si(),    silicon);
  device.set_material(viennashe::materials::si(),    silicon2);
  //
  // Contacts are always conductors/metals (for now)
  device.set_material(viennashe::materials::metal(), gate);
  device.set_material(viennashe::materials::metal(), bulk);


  /**
  Next we set the doping per segment
  In here we use a constant doping for each segment
  and only set a doping on semicondutors
  **/
  std::cout << "* init_device(): Setting doping (per cell) ..." << std::endl;
  device.set_doping_n(Nd, silicon);
  device.set_doping_p(Na, silicon);
  device.set_doping_n(Nd, silicon2);
  device.set_doping_p(Na, silicon2);


  /** Finally, set the contact potentials at the left and the right contacts. **/
  std::cout << "* init_device(): Setting contact potentials (per cell) ..." << std::endl;
  device.set_contact_potential(Vg_init,  gate);
  device.set_contact_potential(0.0,      bulk);

  std::cout << "* init_device(): DONE!" << std::endl;

} // init_device()


/** <h2> The main Simulation Flow</h2>

  With the meshing function generate_mos1d_mesh() and the device
  initialization function init_device() in place, we are ready
  to code up the main application. For simplicity,
  this is directly implemented in the main() routine,
  but a user is free to move this to a separtate function,
  to a class, or whatever other abstraction is appropriate.
  **/
int main()
{
  /** We configure our type of device and use a simple 1D mesh for the MOS-structure **/
  typedef viennashe::device<viennagrid_mesh>    DeviceType;

  std::cout << viennashe::preamble() << std::endl;

  // Create the device object
  DeviceType device;

  /** Before doing anything else, we trigger the device generation: **/
  std::cout << "* main(): Creating mesh ..." << std::endl;

  // Configuration of the 1d MOS
  const double len_gate  =    1e-9;
  const double cs_gate   = 0.01e-9; // Gate cell size
  const double len_oxide =    1e-9;
  const double cs_ox     = 0.05e-9; // Oxide cell size
  const double len_bulk  = 100.e-9;
  const double cs_bulk   =    1e-9; // Bulk cell size

  double Vg = -0.2; // Gate voltage

  // Generate the mesh
  generate_mos1d_mesh(device, len_gate, cs_gate, len_oxide, cs_ox, len_bulk, cs_bulk);

  /** With the mesh in place, we now initalize the device by setting the materials, dopings and contact potentials **/
  //                        ND     NA
  //init_device(device, Vg,  3e23,  1e8 ); // configure a p-channel MOS


  /** <h3>Drift-Diffusion Simulations</h3>

    We start out with solving the drift-diffusion model including quantum corrections.
    For this we first need to set up a configuration object,
    and use this to create and run the simulator object.
  **/
  std::cout << "* main(): Creating DD simulator..." << std::endl;

  /** <h4>Prepare the Drift-Diffusion Simulator Configuration</h4>

    In the next code snippet we set up the configuration for a
    bipolar drift-diffusion simulation. Although most of the options
    we set below are the default values anyway, we recommend the user
    to always set them manually in order to make the code more self-documenting.
   **/
  viennashe::config dd_cfg;
  dd_cfg.with_holes(true);
  dd_cfg.with_electrons(true);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  // Linear solver configuration: We take a dense solver here ...
  dd_cfg.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);

  // Nonlinear solver: 30 Gummel iterations (Newton is also available)
  dd_cfg.nonlinear_solver().max_iters(30);
  dd_cfg.nonlinear_solver().damping(0.6);

  // We will make use of a first-order quantum correction
  dd_cfg.quantum_correction(true);      // solve density gradient equations
  dd_cfg.with_quantum_correction(true); // apply correction potentials from density gradient


  /** <h4>Create and Run the DD-Simulator</h4>
      With the config in place we can create our simulator object.
      Note that after creating your simulator,
      changes to the config *will not* affect the simulator object anymore.
      The simulator is then started using the member function .run()
    **/
  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;
  // Run the DD simulation
  dd_simulator.run();

  /** <h4>Write DD Simulation Output</h4>
    Although one can access all the computed values directly from sources,
    which is somewhat reasonable for one-dimensional device simulations,
    we write the output to simple data files suitable for Gnuplot (http://www.gnuplot.info/)
  **/
  viennashe::io::write_cell_quantity_for_gnuplot(dd_simulator.dg_pot_p(),         device, "mos1d-dg_dd_dgpot_holes.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(dd_simulator.dg_pot_n(),         device, "mos1d-dg_dd_dgpot_electrons.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(dd_simulator.potential(),        device, "mos1d-dg_dd_pot.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(dd_simulator.electron_density(), device, "mos1d-dg_dd_electrons.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(dd_simulator.hole_density(),     device, "mos1d-dg_dd_holes.dat");


  /** <h3>Self-Consistent SHE Simulations</h3>

  In here we do essentially the same as for the drift-diffusion case,
  but instead of a continuity equation we use SHE for holes.
  **/

  std::cout << "* main(): Setting up SHE..." << std::endl;

  /** Create a new config for a SHE simulation
  <h4>Set up the SHE simulator configuration</h4>

  First we set up a new configuration object, enable
  electrons and holes, and specify that we want to use
  SHE for electons, but only a simple continuity equation for holes:
  **/
  viennashe::config config;
  // Configure a first-order SHE for holes only ...
  config.with_holes(true);
  config.set_hole_equation(viennashe::EQUATION_SHE);
  // ... and use DD for the electrons
  config.with_electrons(true);
  config.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  config.max_expansion_order(1);

  // Configure the various scattering mechanisms for SHE
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(true);
  config.scattering().impact_ionization().enabled(false);
  config.scattering().electron_electron(false);

  // Configure the linear solver
  config.linear_solver().max_iters(2000);
  config.linear_solver().ilut_drop_tolerance(1e-2);

  // Nonlinear solver configuration: 30 Gummel steps with damping parameter of 0.5
  config.nonlinear_solver().max_iters(30);
  config.nonlinear_solver().damping(0.5);

  // We require a minimum kinetic energy range of 0.5 eV at each grid node
  config.min_kinetic_energy_range(viennashe::physics::constants::q * 0.5);

  // Set the energy spacing to 31 meV
  config.energy_spacing(viennashe::physics::constants::q * 0.031);

  // Like for DD use density gradient as a first-order quantum correction model
  config.quantum_correction(true);      // calculate a correction potential
  config.with_quantum_correction(true); // use the correction potential for the transport model

  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function run()
   **/
  std::cout << "* main(): Computing SHE..." << std::endl;

  viennashe::simulator<DeviceType> she_simulator(device, config);

  she_simulator.set_initial_guess(viennashe::quantity::potential(), dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(), dd_simulator.hole_density());

  // Run the simulator
  she_simulator.run();

  /** <h4>Write SHE Simulation Output</h4>

  Since the simulation output in (x, H)-space is always at least two-dimensional,
  we first write the results to a VTK file for inspection using e.g. ParaView.
  **/
  std::cout << "* main(): Writing SHE result..." << std::endl;

  /*if (config.get_electron_equation() == viennashe::EQUATION_SHE )
  {
    viennashe::io::she_vtk_writer<DeviceType>()(device,
                                                she_simulator.config(),
                                                she_simulator.quantities().electron_distribution_function(),
                                                "mos1d-dg_edf_n");
  }
  if (config.get_hole_equation() == viennashe::EQUATION_SHE)
  {
    viennashe::io::she_vtk_writer<DeviceType>()(device,
                                                she_simulator.config(),
                                                she_simulator.quantities().hole_distribution_function(),
                                                "mos1d-dg_edf_p");
  }*/

  /**   The macroscopic quantities are written to simple data files suitable for Gnuplot (http://www.gnuplot.info/) **/
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.dg_pot_p(),         device, "mos1d-dg_she_dgpot_holes.dat"     );
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.dg_pot_n(),         device, "mos1d-dg_she_dgpot_electrons.dat" );
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.potential(),        device, "mos1d-dg_she_pot.dat"       );
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.electron_density(), device, "mos1d-dg_she_electrons.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.hole_density(),     device, "mos1d-dg_she_holes.dat");

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e12 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}
