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

// ViennaGrid mesh configurations and centroid() algorithm:
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/algorithm/centroid.hpp"

/** \example np-diode-bipolar.cpp

Two-dimensional bipolar simulation of a np-diode on a structured rectangular grid including generation and recombination are carried out in this example.
An np-diode is the simplest truly bipolar structure and often used for one-dimensional device simulations.
In this example we conduct a two-dimensional device simulation in y-direction for demonstration purposes.
The schematic of the diode is as follows: \code

  -----------------------------------
  |       |      |   |      |       |
  | Metal |  n+  | i |  p+  | Metal |
  |       |      |   |      |       |
  -----------------------------------
  #   0      1     2    3       4        (segment IDs)

\endcode

The intrinsic layer is only one cell thick and only used to improve numerical stability.
The total device length is 225 nanometer with a cell size of 5 nanometer.

<h2>First Step: Initialize the Device</h2>

Before dealing with the actual mesh generation, let's first assume
that we already have the mesh and only need to associate the various
device segments (aka. submeshes) with material parameters, contact
voltages, etc. For simplicity, we collect this initialization
in a separate function:

**/
template <typename DeviceType>
void init_device(DeviceType & device)
{
  typedef typename DeviceType::segment_type          SegmentType;

  /** Define concenience references to the segments: **/
  SegmentType const & contact_left  = device.segment(0);
  SegmentType const & n_left        = device.segment(1);
  SegmentType const & i_center      = device.segment(2);
  SegmentType const & p_right       = device.segment(3);
  SegmentType const & contact_right = device.segment(4);

  /** Set the materials per segment according to the schematic above: **/
  device.set_material(viennashe::materials::metal(), contact_left);
  device.set_material(viennashe::materials::si(),          n_left);
  device.set_material(viennashe::materials::si(),          i_center);
  device.set_material(viennashe::materials::si(),          p_right);
  device.set_material(viennashe::materials::metal(), contact_right);

  /**
  Next we set the doping per segment
  In here we use a constant doping for each segment
  and only set a doping on semicondutors
  **/
  device.set_doping_n(1e22, n_left);
  device.set_doping_p(1e10, n_left);

  device.set_doping_n(1e16, i_center);
  device.set_doping_p(1e16, i_center);

  device.set_doping_n(1e10, p_right);
  device.set_doping_p(1e22, p_right);


  /**
  To also demonstrate the ability of dealing with traps defined over the deivce,
  we add two trap levels (cf. Shockley-Reed-Hall theory) for all silicon segments.
  Note that all values need to be provided in SI units.
  The trap energy is specified relative to the center of the band gap.
  **/
  viennashe::trap_level trap;
  trap.collision_cross_section(3.2e-20);
  trap.density(1e21);
  trap.energy(0.3 * viennashe::physics::constants::q);  //  0.3 eV above center of band gap
  device.add_trap_level(trap);

  // second trap level with different energy, but same density and collision cross section:
  trap.energy(-0.3 * viennashe::physics::constants::q); // -0.3 eV below center of band gap
  device.add_trap_level(trap);


  /** Finally, set the contact potentials at the left and the right contacts. **/
  device.set_contact_potential( 0.0, contact_left);
  device.set_contact_potential(-0.2, contact_right);
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
  /** First we define the device type including the topology to use.
      Here we select a ViennaGrid mesh consisting of quadrilaterals.
      See \ref manual-page-api or the ViennaGrid manual for other mesh types.
   **/
  typedef viennagrid::quadrilateral_2d_mesh                     MeshType;
  typedef viennashe::device<MeshType>                           DeviceType;

  std::cout << viennashe::preamble() << std::endl;

  /** With the device type available, we can directly instantiate an empty device object **/
  std::cout << "* main(): Creating device..." << std::endl;
  DeviceType device;

  /** <h3> Generate the Mesh using the built-in Mesh Generator </h3>
  To generate a suitable mesh for our nin-diode, we use the built-in 1d/2d mesh generator in ViennaSHE.
  In this example we populate the mesh generator configuration accordingly and then pass it to the device,
  which then generates the mesh automatically.
  **/

  viennashe::util::device_generation_config generator_params;

  /** The segments are added one after another via the member function .add_segment().
      For two-dimensional meshes, pass the following six arguments:
        - the x-coordinate of the left-most point
        - the length of the segment in x-direction
        - the number of mesh points to use in x-direction for the segment
        - the y-coordinate of the lowest point
        - the length of the segment in y-direction
        - the number of mesh points to use in y-direction for the segment

  This way the segments are specified one after another:
  **/
  generator_params.add_segment(0.0,    1.0e-7, 2,   //start at x=0,     length 100nm,  2 points
                               0.0,    1.0e-8, 3);  //start at y=0,     length  10nm,  3 points
  generator_params.add_segment(0.0,    1.0e-7, 2,   //start at x=0,     length 100nm,  2 points
                               1.0e-8, 1.0e-7, 21); //start at y= 20nm, length 100nm, 21 points
  generator_params.add_segment(0.0,    1.0e-7, 2,   //start at x=0,     length 100nm,  2 points
                               1.1e-7, 5.0e-9, 2);  //start at y=220nm, length   5nm,  2 points
  generator_params.add_segment(0.0,    1.0e-7, 2,   //start at x=0,     length 100nm,  2 points
                               1.15e-7,1.0e-7, 21); //start at y=230nm, length 100nm, 21 points
  generator_params.add_segment(0.0,    1.0e-7, 2,   //start at x=0,     length 100nm,  2 points
                               2.15e-7,1.0e-8, 3);  //start at y=430,   length  10nm,  3 points

  /** Finally, the mesh generation is triggered by passing the mesh generation to the member function .generate_mesh(): **/
  device.generate_mesh(generator_params);

  /** With the mesh in place, we now initalize the device by setting the materials, dopings and contact potentials.
      This is achieved by the function init_device() we defined above:
  **/
  std::cout << "* main(): Initializing device..." << std::endl;
  init_device(device);

  /** <h3>Drift-Diffusion Simulations</h3>

    We start out with solving the classical drift-diffusion model.
    For this we first set up a configuration object,
    and use this to create and run the simulator object.
  **/
  std::cout << "* main(): Creating DD simulator..." << std::endl;

  /** <h4>Prepare the Drift-Diffusion Simulator Configuration</h4>

    In the next code snippet we set up the configuration for a
    bipolar drift-diffusion simulation. Although most of the options
    we set below are the default values anyway, we recommend the user
    to always set them manually in order to make the code more self-documenting.
   **/
  viennashe::config dd_cfg; // Per default the config is set to DD with Gummel iterations

  // Nonlinear solver: Use up to 100 Gummel iterations with a rather strong damping of 0.2
  dd_cfg.nonlinear_solver().max_iters(100);
  dd_cfg.nonlinear_solver().damping(0.2);

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
    which is somewhat reasonable only for one-dimensional device simulations,
    we write all simulation results to a VTK file, where it can be inspected with e.g. ParaView:
  **/
  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "np-diode_dd_quan");


  /** <h3>Self-Consistent SHE Simulations including Traps</h3>

  In here we do essentially the same as for the drift-diffusion case,
  but instead of a continuity equation we use SHE for both carrier types.
  Also, we enable the consideration of traps explicitly.
  **/
  std::cout << "* main(): Setting up SHE..." << std::endl;

  /** <h4>Set up the SHE Simulator Configuration</h4>

  First we set up a new configuration object, enable
  electrons and holes, and specify that we want to use
  SHE for electrons and holes, and enable traps
  **/
  viennashe::config config;
  // Configure a first-order bipolar SHE
  config.with_electrons(true);
  config.set_electron_equation(viennashe::EQUATION_SHE);

  config.with_holes(true);
  config.set_hole_equation(viennashe::EQUATION_SHE);

  config.with_traps(true);
  config.with_trap_selfconsistency(true);

  // Linear solver: Set up to 1000 iterations
  config.linear_solver().max_iters(1000);

  // Nonlinear solver: Set up to 20 Gummel iterations with a moderate damping of 0.3
  config.nonlinear_solver().max_iters(20);
  config.nonlinear_solver().damping(0.3);

  // SHE: Set first-order expansion and use an energy spacing of 12.5 meV
  config.max_expansion_order(1);
  config.energy_spacing(viennashe::physics::constants::q / 80.0);


  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function .run()
   **/
  std::cout << "* main(): Computing first-order SHE..." << std::endl;

  viennashe::simulator<DeviceType> she_simulator(device, config);

  // Use the DD solution as initial guess
  she_simulator.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());

  she_simulator.run();

  /** <h4>Write SHE Simulation Output</h4>

  The simulation output in (x, H)-space is three-dimensional here,
  so we first write the results to a VTK file for inspection using e.g. ParaView.
  **/
  std::cout << "* main(): Writing SHE result..." << std::endl;
  viennashe::io::she_vtk_writer<DeviceType>()(device,
                                              she_simulator.config(),
                                              she_simulator.quantities().electron_distribution_function(),
                                              "np-diode-she");

  /** Moreover, we write macroscopic quantities, in particular the average trap occupancy, to a another VTK file:
  **/
  viennashe::io::write_quantities_to_VTK_file(she_simulator, "np-diode_she_quan");

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e12 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}

