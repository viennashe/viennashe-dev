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


/** \example nin-diode-1d.cpp

  Here we consider the simulation of a purely 1d nin-diode.
  This example also illustrates the structured 1d device generator.

  A schematic of the nin-diode, which is one of the simplest structures and often used for one-dimensional device simulations, is as follows:

  \code

    -------------------------------------
    |       |      |     |      |       |
    | Metal |  n+  |  n  |  n+  | Metal |
    |       |      |     |      |       |
    -------------------------------------
    #   0       1     2     3       4        (segment IDs)

  \endcode

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
  SegmentType contact_left  = device.segment(0);
  SegmentType n_left        = device.segment(1);
  SegmentType i_center      = device.segment(2);
  SegmentType n_right       = device.segment(3);
  SegmentType contact_right = device.segment(4);

  /** Set the materials per segment according to the schematic above: **/
  device.set_material(viennashe::materials::si(),    n_left);
  device.set_material(viennashe::materials::si(),    i_center);
  device.set_material(viennashe::materials::si(),    n_right);

  // Contacts are always conductors/metals (for now)
  device.set_material(viennashe::materials::metal(), contact_left);
  device.set_material(viennashe::materials::metal(), contact_right);

  /**
  Next we set the doping per segment
  In here we use a constant doping for each segment
  and only set a doping on semicondutors
  **/
  std::cout << "* init_device(): Setting doping..." << std::endl;
  device.set_doping_n(1e25, n_left);
  device.set_doping_p(1e6,  n_left);

  device.set_doping_n(1e20, i_center);
  device.set_doping_p(1e11, i_center);

  device.set_doping_n(1e25, n_right);
  device.set_doping_p(1e6,  n_right);

  /** Finally, set the contact potentials at the left and the right contacts. **/
  device.set_contact_potential(0.0, contact_left);
  device.set_contact_potential(0.5, contact_right);
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
  typedef viennashe::device            DeviceType;

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

  /** The segments are then added one after another via the member function .add_segment().
      For one dimensional meshes, pass the following three arguments:
        - the x-coordinate of the left-most point
        - the length of the segment
        - the number of mesh points to use for the segment

  This way the segments are specified one after another:
  **/
  // Left contact
  generator_params.add_segment(0.0,  1e-6, 10);  //start at x=0, length 1e-6, 10 points
  // n_left
  generator_params.add_segment(1e-6, 1e-6, 10);
  // i_center
  generator_params.add_segment(2e-6, 1e-6, 10);
  // n_right
  generator_params.add_segment(3e-6, 1e-6, 10);
  // Right contact
  generator_params.add_segment(4e-6, 1e-6, 10);

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
  viennashe::config dd_cfg;

  // enable electrons and holes, use the continuity equation for each of them:
  dd_cfg.with_holes(true);
  dd_cfg.with_electrons(true);

  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);

  // The default nonlinear solver is a Gummel iteration. Here we use Newton's method with at most 50 iterations:
  dd_cfg.nonlinear_solver().set(viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);
  dd_cfg.nonlinear_solver().max_iters(50);

  /** <h4>Create and Run the DD-Simulator</h4>
      With the config in place we can create our simulator object.
      Note that after creating your simulator,
      changes to the config *will not* affect the simulator object anymore.
      The simulator is then started using the member function .run()
    **/
  viennashe::simulator dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;
  dd_simulator.run();

  /** <h4>Write DD Simulation Output</h4>
    Although one can access all the computed values directly from sources,
    which is somewhat reasonable for one-dimensional device simulations,
    we write all simulation results to a VTK file, where it can be inspected with e.g. ParaView:
  **/
  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "nin1d_dd_quan");


  /** <h3>Self-Consistent SHE Simulations</h3>

  In here we do essentially the same as for the drift-diffusion case,
  but instead of a continuity equation we use SHE for both carrier types.
  **/
  std::cout << "* main(): Setting up SHE..." << std::endl;

  /** <h4>Set up the SHE Simulator Configuration</h4>

  First we set up a new configuration object, enable
  electrons and holes, and specify that we want to use
  SHE for electons, but only a simple continuity equation for holes:
  **/
  viennashe::config config;

  // Config a first-order SHE for electrons and holes
  config.max_expansion_order(1);
  config.with_electrons(true);
  config.with_holes(true);
  config.set_electron_equation(viennashe::EQUATION_SHE);
  config.set_hole_equation(viennashe::EQUATION_SHE);

  // Scattering mechanisms for SHE: We only consider phonon scattering:
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(false);
  config.scattering().impact_ionization().enabled(false);
  config.scattering().electron_electron(false);

  // nonlinear solver: At most 30 Gummel iterations with a damping of 0.3:
  config.nonlinear_solver().max_iters(30);
  config.nonlinear_solver().damping(0.4);

  // Configure SHE: Minimum energy range is 0.5 eV
  config.min_kinetic_energy_range(viennashe::physics::constants::q * 0.5);

  // energy spacing of 31meV
  config.energy_spacing(viennashe::physics::constants::q * 0.031);


  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function .run()
   **/
  std::cout << "* main(): Computing SHE..." << std::endl;
  viennashe::simulator she_simulator(device, config);

  // Use DD solution as an initial guess
  she_simulator.set_initial_guess(viennashe::quantity::potential(), dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(), dd_simulator.hole_density());

  she_simulator.run();

  /** <h4>Write SHE Simulation Output</h4>

  Since the simulation output in (x, H)-space is always at least two-dimensional,
  we first write the results to a VTK file for inspection using e.g. ParaView.
  **/
  std::cout << "* main(): Writing SHE result..." << std::endl;

  viennashe::io::she_vtk_writer()(device,
                                  she_simulator.config(),
                                  she_simulator.quantities().electron_distribution_function(),
                                  "nin1d_edf");

  /** In addition, we write all macroscopic quantities such as carrier concentrations to a VTK file:
  **/
  viennashe::io::write_quantities_to_VTK_file(she_simulator, "nin1d_she_quan");


  /** To easily plot 1d-graphs showing the computed energy distribution function, we also write the results
      to simple text files, which can be processed by e.g. Gnuplot (http://www.gnuplot.info/) **/
  viennashe::io::gnuplot_edf_writer gedfwriter;

  // Helper type returning the generalized energy distribution function:
  typedef viennashe::she::generalized_edf_wrapper<viennashe::simulator::she_quantity_type> EDFWrapperType;

  // Write the generalized EDF to a plain text file
  gedfwriter(device,
             viennashe::util::any_filter(), // write the result for all cells
             EDFWrapperType(she_simulator.config(), she_simulator.quantities().electron_distribution_function()),
             "nin1d_edf.dat");              // filename

  // Write the electron density to a plain text file
  viennashe::io::gnuplot_writer gpwriter;
  gpwriter(device,
           viennashe::util::any_filter(),   // write the result for all cells
           she_simulator.electron_density(),
           "nin1d_n.dat");                  // filename

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e12 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}

