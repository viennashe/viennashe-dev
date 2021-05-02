/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:      MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include <iostream>
#include <cstdlib>
#include <vector>
#include "iacaMarks.h"

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid mesh configurations:
#include "viennagrid/config/default_configs.hpp"



/** \example mosfet-dg.cpp

  This example is similar to mosfet.cpp, where a two-dimensional MOSFET simulation is carried out.
  This time, however, also quantum corrections as computed from the density gradient model are considered.
  A schematic of the considered device with segment numbers is as follows:
<table>
 <tr>
  <td>
   <table cellspacing="0" cellpadding="0" >
   <tr><th>Segment #</th><th> Segment description</th><th>Notes</th></tr>
   <tr><td> 1</td><td> Gate contact.  </td><td>Potential known.          </td></tr>
   <tr><td> 2</td><td> Source contact.</td><td>Potential known.          </td></tr>
   <tr><td> 3</td><td> Oxide.         </td><td>No boundary conditions, no carriers here.   </td></tr>
   <tr><td> 4</td><td> Drain contact. </td><td>Potential known.          </td></tr>
   <tr><td> 5</td><td> Source region. </td><td>n-doped region.           </td></tr>
   <tr><td> 6</td><td> Body.          </td><td>Intrinsic region.         </td></tr>
   <tr><td> 7</td><td> Drain region.  </td><td>n-doped region.           </td></tr>
   <tr><td> 8</td><td> Bulk contact.  </td><td>Potential known.          </td></tr>
   <tr><td colspan="3">See Netgen geometry description in mosfet.in2d</td></tr>
   </table>
  </td>
  <td>
   ![Schematic of the considered MOSFET](mosfet-schematic.png)
  </td>
 </tr>
</table>

<h2>First Step: Initialize the Device</h2>

Once the device mesh is loaded, we need to initialize the various
device segments (aka. submeshes) with material parameters, contact
voltages, etc. For simplicity, we collect this initialization
in a separate function:
**/
template <typename DeviceType>
void init_device(DeviceType & device)
{
  typedef typename DeviceType::segment_type          SegmentType;

  /** Provide convenience names for the various segments: **/
  SegmentType const & gate_contact   = device.segment(1);
  SegmentType const & source_contact = device.segment(2);
  SegmentType const & gate_oxide     = device.segment(3);
  SegmentType const & drain_contact  = device.segment(4);
  SegmentType const & source         = device.segment(5);
  SegmentType const & drain          = device.segment(6);
  SegmentType const & body           = device.segment(7);
  SegmentType const & body_contact   = device.segment(8);

  /** Now we are ready to set the material for each segment: **/
  std::cout << "* init_device(): Setting materials..." << std::endl;
  device.set_material(viennashe::materials::metal(),   gate_contact);
  device.set_material(viennashe::materials::metal(), source_contact);
  device.set_material(viennashe::materials::metal(),  drain_contact);
  device.set_material(viennashe::materials::metal(),   body_contact);

  device.set_material(viennashe::materials::hfo2(), gate_oxide);

  device.set_material(viennashe::materials::si(),  source);
  device.set_material(viennashe::materials::si(),   drain);
  device.set_material(viennashe::materials::si(),    body);

  /** For all semiconductor cells we also need to specify a doping.
      If the doping is inhomogeneous, one usually wants to set this
      through some automated process (e.g. reading from file).
      For simplicitly we use a doping profile which is constant per segment.
      Note that the doping needs to be provided in SI units, i.e. \f$m^{-3}\f$
  **/
  std::cout << "* init_device(): Setting doping..." << std::endl;
  device.set_doping_n(1e24, source);
  device.set_doping_p(1e8,  source);

  device.set_doping_n(1e24, drain);
  device.set_doping_p(1e8,  drain);

  device.set_doping_n(1e17, body);
  device.set_doping_p(1e15, body);

  /** Finally, we need to provide contact potentials for the device.
      Since we already have dedicated contact segments,
      all we need to do is to set the contact voltages per segment:
  **/
  device.set_contact_potential(0.8, gate_contact);
  device.set_contact_potential(0.0, source_contact);
  device.set_contact_potential(0.1, drain_contact);
  device.set_contact_potential(0.0, body_contact);

}


/** <h2> The main Simulation Flow</h2>

  With the function init_device() in place, we are ready
  to code up the main application. For simplicity,
  this is directly implemented in the main() routine,
  but a user is free to move this to a separate function,
  to a class, or whatever other abstraction is appropriate.
  **/
int main()
{
  /** First we define the device type including the topology to use.
      Here we select a ViennaGrid mesh consisting of triangles.
      See \ref manual-page-api or the ViennaGrid manual for other mesh types.
   **/
  typedef viennashe::device<viennagrid::triangular_2d_mesh>       DeviceType;
  typedef DeviceType::segment_type                                SegmentType;

  std::cout << viennashe::preamble() << std::endl;

  /** <h3>Read and Scale the Mesh</h3>
      Since it is inconvenient to set up a big triangular mesh by hand,
      we load a mesh generated by Netgen. The spatial coordinates
      of the Netgen mesh are in nanometers, while ViennaSHE expects
      SI units (meter). Thus, we scale the mesh by a factor of \f$ 10^{-9} \f$.
  **/

  std::cout << "* main(): Creating and scaling device..." << std::endl;
  DeviceType device;
  device.load_mesh("../examples/data/mosfet840.mesh");
  device.scale(1e-9);

//while(cond)
{
    IACA_START
  /** <h3>Initialize the Device</h3>
    Here we just need to call the initialization routine defined before:
   **/
  std::cout << "* main(): Initializing device..." << std::endl;
  init_device(device);

  /** <h3>Drift-Diffusion Simulations</h3>

    In order to compute a reasonable initial guess of
    the electrostatic potential for SHE,
    we first solve the drift-diffusion model.
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

  // Nonlinear solver parameters: 50 Gummel iterations with moderate damping
  dd_cfg.nonlinear_solver().max_iters(50);
  dd_cfg.nonlinear_solver().damping(0.5);

  /** Quantum Corrections:
  Enable the computation of quantum corrections via the density gradient model: **/
  dd_cfg.quantum_correction(true);
  /** However, Do *not* use the correction potential in the continuity equations and just compute classic carrier concentrations.
      Here we just want to know the correction potential as a kind of error estimate. **/
  dd_cfg.with_quantum_correction(false);

  /** <h4>Create and Run the DD-Simulator</h4>
      With the config in place, we can create our simulator object.
      Note that after creating your simulator object,
      changes to the config *will not* affect the simulator object anymore.
      The simulator is then started using the member function .run()
    **/
  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);
  std::cout << "* main(): Launching DD simulator..." << std::endl;
  dd_simulator.run();

  /** <h4>Write DD Simulation Output</h4>
    Although one can access all the computed values directly from sources,
    for typical meshes this is way too tedious to do by hand.
    Thus, the recommended method for inspecting simulator output
    is by writing the computed values to a VTK file, where
    it can then be inspected by e.g. ParaView.
  **/
  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "mosfet-dg_dd_quan");

  // Write only the density gradient correction to a separate VTK file:
  viennashe::io::write_quantity_to_VTK_file(dd_simulator.dg_pot_n(), device, "mosfet-dg_dgpot");

  /** <h4>Calculate Terminal Currents</h4>
    Since the terminal currents are not directly visible in the VTK files, we compute them directly here.
    To simplify matters, we only output the electron and hole drain currents from the body segment into the drain contact:
  **/
  SegmentType const & drain_contact  = device.segment(4);
  SegmentType const & body           = device.segment(6);

  std::cout << "* main(): Drain electron current Id_e = " << viennashe::get_terminal_current(device, viennashe::ELECTRON_TYPE_ID,
      dd_simulator.potential(), dd_simulator.electron_density(),
      viennashe::models::create_constant_mobility_model(device, 0.1430),
      body, drain_contact ) * 1e-6 << std::endl;
  std::cout << "* main(): Drain hole current Id_h = " << viennashe::get_terminal_current(device, viennashe::HOLE_TYPE_ID,
      dd_simulator.potential(), dd_simulator.hole_density(),
      viennashe::models::create_constant_mobility_model(device, 0.0460),
      body, drain_contact ) * 1e-6 << std::endl;

  /** <h3>Self-Consistent SHE Simulations</h3>
    To run self-consistent SHE simulations, we basically proceed as for the drift-diffusion case above,
    but have to explicitly select the SHE equations.

    Similar to the case of the drift-diffusion simulation above,
    we first need to set up the configuration.

    <h4>Prepare the SHE simulator configuration</h4>

    First we set up a new configuration object, enable
    electrons and holes, and specify that we want to use
    SHE for electons, but only a simple continuity equation for holes:
  **/
  std::cout << "* main(): Setting up first-order SHE (semi-self-consistent using 40 Gummel iterations)..." << std::endl;
  viennashe::config config;

  // Set the expansion order to 1 for a first-order SHE
  config.max_expansion_order(1);

  // Use both carrier types
  config.with_electrons(true); config.with_holes(true);

  // Configure equations: SHE for electrons, continuity for holes:
  config.set_electron_equation(viennashe::EQUATION_SHE); // SHE for electrons
  config.set_hole_equation(viennashe::EQUATION_CONTINUITY); // DD for holes

  // Energy spacing of 31 meV
  config.energy_spacing(31.0 * viennashe::physics::constants::q / 1000.0);

  // Set the scattering mechanisms: We use phonon scattering, but no ionized impurity scattering
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(false);

  // The linear solver should run for at most 2000 iterations:
  config.linear_solver().max_iters(2000);

  // The nonlinear solver should use up to 40 Gummel iterations with a damping parameter of 0.4:
  config.nonlinear_solver().max_iters(40);
  config.nonlinear_solver().damping(0.4);

  // Configure the first-order quantum correction
  config.quantum_correction(true);       // Solve density gradient ...
  config.with_quantum_correction(false); // ... but don't use it to correct the potential (cf. DD solver above)

  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function run()
   **/
  std::cout << "* main(): Computing first-order SHE..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);

  // Set the previous DD solution as an initial guess
  she_simulator.set_initial_guess(viennashe::quantity::potential(), dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(), dd_simulator.hole_density());
  she_simulator.set_initial_guess(viennashe::quantity::density_gradient_electron_correction(), dd_simulator.dg_pot_n());
  she_simulator.set_initial_guess(viennashe::quantity::density_gradient_hole_correction(), dd_simulator.dg_pot_p());

  // Run the simulation
  she_simulator.run();

  /** <h4>Write SHE Simulation Output</h4>

  With a spatially two-dimensional mesh, the result in (x, H)-space is three-dimensional.
  The solutions computed in this augmented space are written to a VTK file for inspection using e.g. ParaView:
  **/
  std::cout << "* main(): Writing energy distribution function from first-order SHE result..." << std::endl;
  viennashe::io::she_vtk_writer<DeviceType>()(device,
                                              she_simulator.config(),
                                              she_simulator.quantities().electron_distribution_function(),
                                              "mosfet-dg_she_edf");

  /** Here we also write the potential and electron density to separate VTK files: **/
  viennashe::io::write_quantity_to_VTK_file(she_simulator.potential(), device, "mosfet-dg_she_potential");
  viennashe::io::write_quantity_to_VTK_file(she_simulator.electron_density(), device, "mosfet-dg_she_electrons");

  /** Write all macroscopic result quantities (carrier concentrations, density gradient corrections, etc.) to a single VTK file:
  **/
  viennashe::io::write_quantities_to_VTK_file(she_simulator, "mosfet-dg_she_quan");

  /** <h4>Calculate Terminal Currents</h4>
    Since the terminal currents are not directly visible in the VTK files, we compute them directly here.
    To simplify matters, we only compute the electron current from the body segment into the drain contact based on the solution of the SHE equations:
  **/
  std::cout << "* main(): Drain electron current Id_e = " << viennashe::get_terminal_current(
      device, config, she_simulator.quantities().electron_distribution_function(), body, drain_contact ) * 1e-6
      << std::endl;

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e12 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;
}
IACA_END

  return EXIT_SUCCESS;
}

