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
#include "viennagrid/config/default_configs.hpp"


/** \example resistor.cpp

  The simulation of a resistor is probably the simplest device one could think of.
  We use it to run a classical drift-diffusion simulation, followed by a first-order SHE simulation,
  and finally we run a SHE simulation including electron-electron scattering.

  The device schematic is as follows:
  \code
    -------------------------
    | |                   | |
    |M|      Silicon      |M|
    | |                   | |
    -------------------------
  \endcode
  where the first and the last cell are used as contacts.

<h2>First Step: Initialize the Device</h2>

Before dealing with the actual mesh generation, let's first assume
that we already have the mesh and only need to associate the various
device segments (aka. submeshes) with material parameters, contact
voltages, etc. For simplicity, we collect this initialization
in a separate function:
**/
template <typename DeviceType>
void init_device(DeviceType & device, double vcc)
{
  // ViennaGrid types for mesh handling:
  typedef typename DeviceType::mesh_type           MeshType;

  typedef typename viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;

  // Container of all cells in the mesh:
  CellContainer cells(device.mesh());

  /** Set the materials per segment according to the schematic above: **/
  device.set_material(viennashe::materials::si()); // Silicon everywhere

  // overwrite the material of the two cells on the left and right:
  device.set_material(viennashe::materials::metal(), cells[0]);
  device.set_material(viennashe::materials::metal(), cells[cells.size()-1]);

  /**
  Next we set the doping. Since the device is homogeneous, we can just specify it throughout the device:
  Note that the doping in metal regions is ignored.
  Also, keep in mind that the doping must be provided in SI units, i.e. per cubic-meter.
  **/
  std::cout << "* init_device(): Setting doping..." << std::endl;
  device.set_doping_n(1e28);
  device.set_doping_p(1e4);

  /** Set the contact potentials: **/
  device.set_contact_potential(0.0, cells[0]);
  device.set_contact_potential(vcc, cells[cells.size()-1]);
}


/** <h2> The main Simulation Flow</h2>

  With the device initialization function init_device() in place,
  we are ready to code up the main application. For simplicity,
  this is directly implemented in the main() routine,
  but a user is free to move this to a separate function,
  to a class, or whatever other abstraction is appropriate.
  **/
int main()
{
  /** We configure our device type using a simple 1D mesh consisting of lines.
      In one dimension there is no other element type available, while in
      two and three spatial dimensions one can choose from triangles, quadrilaterals, etc.
  **/
  typedef viennagrid::line_1d_mesh       MeshType;
  typedef viennashe::device<MeshType>    DeviceType;

  std::cout << viennashe::preamble() << std::endl;

  /** With the device type available, we can directly instantiate an empty device object **/
  DeviceType device;


  /** <h3> Generate the Mesh using the built-in Mesh Generator </h3>
  To generate a suitable mesh for our nin-diode, we use the built-in 1d/2d mesh generator in ViennaSHE.
  In this example we populate the mesh generator configuration accordingly and then pass it to the device,
  which then generates the mesh automatically.
  **/
  viennashe::util::device_generation_config generator_params;

  /** Here we only need one segment, for which we pass the following three arguments:
        - the x-coordinate of the left-most point
        - the length of the segment
        - the number of mesh points to use for the segment

  **/
  generator_params.add_segment(0.0,  1e-6, 11);  //start at x=0, length 1um, 11 points

  /** Finally, the mesh generation is triggered by passing the mesh generation to the member function .generate_mesh(): **/
  device.generate_mesh(generator_params);


  /** With the mesh in place, we now initalize the device by setting the materials, dopings and contact potentials.
      This is achieved by the function init_device() we defined above:
  **/
  std::cout << "* main(): Initializing device..." << std::endl;
  init_device(device, 1.0);


  /** <h3>Drift-Diffusion Simulations</h3>

    We start out with solving the classical drift-diffusion model.
    For this we first set up a configuration object,
    and use it to create and run the simulator object.
  **/
  std::cout << "* main(): Creating DD simulator..." << std::endl;

  /** <h4>Prepare the Drift-Diffusion Simulator Configuration</h4>

    In the next code snippet we set up the configuration for a
    bipolar drift-diffusion simulation. Although most of the options
    we set below are the default values anyway, we recommend the user
    to always set them manually in order to make the code more self-documenting.
   **/
  viennashe::config dd_cfg;

  // Linear solver: Use a dense linear solver, because the mesh is small (no preconditioning troubles):
  dd_cfg.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);

  // Nonlinear solver: 20 Gummel iterations with a damping parameter of 0.25:
  dd_cfg.nonlinear_solver().max_iters(20);
  dd_cfg.nonlinear_solver().damping(0.25);

  /** <h4>Create and Run the DD-Simulator</h4>
      With the config in place we can create our simulator object.
      Note that after creating your simulator,
      changes to the config *will not* affect the simulator object anymore.
      The simulator is then started using the member function .run()
    **/
  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;
  dd_simulator.run();


  /** <h4>Write DD Simulation Output</h4>
    Although one can access all the computed values directly from sources,
    which is somewhat reasonable for one-dimensional device simulations,
    we write all simulation results to a VTK file, where it can be inspected with e.g. ParaView:
  **/
  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "resistor_dd_quan");

  /** <h3>Fixed-Field SHE Simulations</h3>

  In here we do essentially the same as for the drift-diffusion case,
  but instead of a continuity equation we use SHE for electrons.
  Since we want the potential to be linear along the device without any spurious
  effects at contacts, the simulation is not self-consistent.
  **/
  std::cout << "* main(): Setting up SHE..." << std::endl;

  /** <h4>Set up the SHE Simulator Configuration</h4>

  First we set up a new configuration object, enable
  electrons and holes, and specify that we want to use
  SHE for electons, but only a simple continuity equation for holes:
  **/
  viennashe::config config;

  config.with_electrons(true);
  config.set_electron_equation(viennashe::EQUATION_SHE);

  config.with_holes(true);
  config.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  // Nonlinear solver: Only one Gummel iteration with very high damping (to avoid changes to the potential)
  config.nonlinear_solver().max_iters(1);
  config.nonlinear_solver().damping(0.001);

  // SHE: First-order simulation using an energy spacing of 15.5meV and a kinetic energy range of at least 1eV:
  config.max_expansion_order(1);
  config.energy_spacing(0.0155 * viennashe::physics::constants::q);
  config.min_kinetic_energy_range(1.0 * viennashe::physics::constants::q);

  // configure scattering mechanisms. We use phonon scattering only for now.
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(false);
  config.scattering().impact_ionization().enabled(false);
  config.scattering().electron_electron(false);


  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function .run()
   **/
  std::cout << "* main(): Computing SHE (electrons only) ..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);

  she_simulator.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());

  she_simulator.run();

  std::cout << "* main(): Writing SHE result..." << std::endl;

  /** <h4>Write SHE Simulation Output</h4>

  Write the computed distribution function to a plain text file,
  from where it can be processed using e.g. Gnuplot (http://www.gnuplot.info/)
  **/
  viennashe::io::gnuplot_edf_writer edfwriter;
  edfwriter(device, viennashe::util::any_filter(), she_simulator.edf(viennashe::ELECTRON_TYPE_ID), "resistor_edf.dat");

  /** Also write the current density and the electron density: **/
  typedef viennashe::simulator<DeviceType>::she_quantity_type she_quan_type;
  viennashe::she::current_density_wrapper<DeviceType, she_quan_type> Jn(device, config,
                she_simulator.quantities().electron_distribution_function());

  viennashe::io::write_cell_quantity_for_gnuplot(Jn, device, "resistor_Jn.dat");
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator.electron_density(), device, "resistor_n.dat");

  /** Write all computed quantities to VTK for inspection with e.g. ParaView: **/
  viennashe::io::write_quantities_to_VTK_file(she_simulator, "resistor_she_quan");


  /** <h3>SHE Simulations including Electron-Electron Scattering</h3>

  So far the simulations only included the interaction of carriers with the crystal lattice.
  If carriers gain high energies, then electron-electron scattering has a significant
  influence on the high-energy tail of the distribution function.
  Here we study the impact of electron-electron scattering on the resistor.

  <h4>Set up the SHE Simulator Configuration</h4>
  For simplicity we simply reuse the previous configuration.
  The only difference is that we allow for 10 Gummel iterations
  (because electron-electron scattering is a nonlinear effect)
  and enable the process explicitly:
  **/
  std::cout << "* main(): Computing SHE (electrons only) with EE scattering ..." << std::endl;
  config.nonlinear_solver().max_iters(10);
  config.scattering().electron_electron(true);

  /** <h4>Create and Run the SHE-Simulator</h4>
    Just as before, the simulator object is instantiated by passing the device and the configuration to the constructor.
    The potential as well as carrier concentrations are taken from the previous first-order SHE simulation:
   **/
  viennashe::simulator<DeviceType> she_simulator_ee(device, config);

  she_simulator_ee.set_initial_guess(viennashe::quantity::potential(),        she_simulator.potential());
  she_simulator_ee.set_initial_guess(viennashe::quantity::electron_density(), she_simulator.electron_density());
  she_simulator_ee.set_initial_guess(viennashe::quantity::hole_density(),     she_simulator.hole_density());

  she_simulator_ee.run();

  /** <h4>Write SHE Simulation Output</h4>

  As before, we write the computed distribution function to a plain text file,
  from where it can be processed using e.g. Gnuplot (http://www.gnuplot.info/)
  **/
  std::cout << "* main(): Writing SHE result..." << std::endl;

  edfwriter(device, viennashe::util::any_filter(), she_simulator_ee.edf(viennashe::ELECTRON_TYPE_ID), "resistor_edf_ee.dat");

  /** Let's also write the electron density and the average carrier energy to a plain text file:
  **/
  viennashe::io::write_cell_quantity_for_gnuplot(she_simulator_ee.electron_density(), device, "resistor_n_ee.dat");

  viennashe::she::carrier_energy_wrapper<she_quan_type> avge (config, she_simulator_ee.quantities().electron_distribution_function());
  viennashe::io::write_cell_quantity_for_gnuplot(avge, device, "resistor_avg_ee.dat");

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e12 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}

