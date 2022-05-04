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

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include <iostream>
#include <cstdlib>
#include <vector>

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default configurations:
#include "viennagrid/config/default_configs.hpp"


/** \example half-trigate.cpp

 In this example we simulate a 3d FinFET (trigate) transistor, which is slit along the plane of symmetry.
 The device mesh is partitioned as follows:
<table>
 <tr>
  <td>
   <table cellspacing="0" cellpadding="0" >
   <tr><th>Segment #</th><th> Segment description</th><th>Notes</th></tr>
   <tr><td>1</td><td> Source.         </td><td>Potential known at contact.   </td></tr>
   <tr><td>2</td><td> Channel.        </td><td>No boundary conditions.       </td></tr>
   <tr><td>3</td><td> Drain.          </td><td>Potential known at contact.   </td></tr>
   <tr><td>4</td><td> Oxide.          </td><td>No carriers here.             </td></tr>
   <tr><td>5</td><td> Gate.           </td><td>Potential known. Thermal bath.</td></tr>
   <tr><td>6</td><td> Body.           </td><td>Potential known at bottom.    </td></tr>
   <tr><td>7</td><td> Source Contact. </td><td>Potential known. Thermal bath.</td></tr>
   <tr><td>8</td><td> Drain Contact.  </td><td>Potential known. Thermal bath.</td></tr>
   <tr><td>9</td><td> Body Contact.   </td><td>Potential known. Thermal bath.</td></tr>
   <tr><td colspan="3">See Netgen geometry description in trigate.geo</td></tr>
   </table>
  </td>
  <td>
 ![Half of a FinFET/trigate device (symmetric)](half-trigate.png)
  </td>
  <td>
  ![Average Carrier Energy (eV)](half-trigate-energy.png)
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
  SegmentType const &  source_segment = device.segmentation()[1];
  SegmentType const & channel_segment = device.segmentation()[2];
  SegmentType const &   drain_segment = device.segmentation()[3];
  SegmentType const &   oxide_segment = device.segmentation()[4];
  SegmentType const &    gate_segment = device.segmentation()[5];
  SegmentType const &    body_segment = device.segmentation()[6];

  SegmentType const &  source_contact_segment = device.segmentation()[7];
  SegmentType const &   drain_contact_segment = device.segmentation()[8];
  SegmentType const &    body_contact_segment = device.segmentation()[9];

  /** Now we are ready to set the material for each segment: **/

  std::cout << "* init_device(): Setting materials..." << std::endl;
  device.set_material(viennashe::materials::metal(), gate_segment);
  device.set_material(viennashe::materials::hfo2(), oxide_segment);
  device.set_material(viennashe::materials::si(),  source_segment);
  device.set_material(viennashe::materials::si(), channel_segment);
  device.set_material(viennashe::materials::si(),   drain_segment);
  device.set_material(viennashe::materials::si(),    body_segment);

  /** Note that contacts are always conductors/metals (for now) **/
  device.set_material(viennashe::materials::metal(), source_contact_segment);
  device.set_material(viennashe::materials::metal(),  drain_contact_segment);
  device.set_material(viennashe::materials::metal(),   body_contact_segment);

  /** For all semiconductor cells we also need to specify a doping.
      If the doping is inhomogeneous, one usually wants to set this
      through some automated process (e.g. reading from file).
      Here we use a doping profile which is constant per segment.
      Note that the doping needs to be provided in SI units, i.e. \f$m^{-3}\f$
  **/
  std::cout << "* init_device(): Setting doping..." << std::endl;

  device.set_doping_n(1e18, channel_segment);
  device.set_doping_p(1e14, channel_segment);

  device.set_doping_n(1e26, source_segment);
  device.set_doping_p(1e6,  source_segment);

  device.set_doping_n(1e26, drain_segment);
  device.set_doping_p(1e6,  drain_segment);

  device.set_doping_n(1e18, body_segment);
  device.set_doping_p(1e14, body_segment);

  /** Finally, we need to provide contact potentials for the device.
      Since we already have dedicated contact segments,
      all we need to do is to set the contact voltages per segment:
  **/
  std::cout << "* init_device(): Setting contacts..." << std::endl;

  double gnd   = 0.0;
  double vcc   = 0.3;
  double vgate = 0.8;

  device.set_contact_potential(vgate, gate_segment);
  device.set_contact_potential(gnd,   source_contact_segment);
  device.set_contact_potential(vcc,   drain_contact_segment);
  device.set_contact_potential(gnd,   body_contact_segment);
}

/** <h2> The main Simulation Flow</h2>

  With the function init_device() in place, we are ready
  to code up the main application. For simplicity,
  this is directly implemented in the main() routine,
  but a user is free to move this to a separtate function,
  to a class, or whatever other abstraction is appropriate.
  **/
int main()
{
  /** First we define the device type including the topology to use.
      Here we select a ViennaGrid mesh consisting of tetrahedra.
      See \ref manual-page-api or the ViennaGrid manual for other mesh types.
   **/
  typedef viennashe::device<viennagrid::tetrahedral_3d_mesh> DeviceType;

  std::cout << viennashe::preamble() << std::endl;

  /** <h3>Read and Scale the Mesh</h3>
      Since it is inconvenient to set up a tetrahedral mesh by hand,
      we load a mesh generated by Netgen. The spatial coordinates
      of the Netgen mesh are in nanometers, while ViennaSHE expects
      SI units (meter). Thus, we scale the mesh by a factor of \f$ 10^{-9} \f$.
  **/
  std::cout << "* main(): Creating and scaling device..." << std::endl;
  DeviceType device;
  try
  {
    device.load_mesh("../examples/data/half-trigate57656.mesh");
  }
  catch (std::runtime_error const & e)
  {
    std::cerr << "-----------------------------------------------------------" << std::endl;
    std::cerr << "--- NOTE: Please download the mesh file from:" << std::endl;
    std::cerr << "--- http://viennashe.sourceforge.net/half-trigate57656.mesh" << std::endl;
    std::cerr << "--- and place it in folder ../examples/data/" << std::endl;
    std::cerr << "--- (this way the repository remains light-weight)" << std::endl;
    std::cerr << "-----------------------------------------------------------" << std::endl;
    throw e;
  }

  device.scale(1e-9);


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

  // Create the configuration object
  viennashe::config dd_cfg;

  // Specify that the drift-diffusion system for electrons and holes,
  // both using the continuity equations, should be solved:
  dd_cfg.with_electrons(true);
  dd_cfg.with_holes(true);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  /** For the non-linear solver we use the Gummel method.
     To enable a Newton-Raphson solver one could use use:
  **/
  //   dd_cfg.nonlinear_solver().set(viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

  /** In this example we stick with the default Gummel method.
      Typically, the Gummel iteration without a damping will diverge,
      hence we set a damping parameter. Values closer to zero
      result in a stronger damping, while a damping value of unity
      implies no damping at all. Values between 0.1 and 0.5 work in
      many cases. On the other hand, lower values for the damping parameter
      require higher iteration counts, hence the two values should be adjusted
      concurrently. **/
  dd_cfg.nonlinear_solver().damping(0.4);  // The damping factor
  dd_cfg.nonlinear_solver().max_iters(50); // Maximum of nonlinear iterations

  /** <h4>Create and Run the DD-Simulator</h4>
      With the config in place we can create our simulator object.
      Note that after creating your simulator,
      changes to the config *will not* affect the simulator object anymore.
      The simulator is then started using the member function .run()
    **/
  std::cout << "* main(): Creating and launching DD simulator..." << std::endl;
  viennashe::simulator<DeviceType> dd_simulator(device,   // the device created before
                                                dd_cfg);  // the configuration
  dd_simulator.run();

  /** <h4>Write DD Simulation Output</h4>
    Although one can access all the computed values directly from sources,
    for typical meshes this is way too tedious to do by hand.
    Thus, the recommended method for inspecting simulator output
    is by writing the computed values to a VTK file, where
    it can then be inspected by e.g. ParaView.
  **/
  viennashe::io::write_quantities_to_VTK_file(dd_simulator,            // simulator object
                                              "half-trigate_dd_quan"); // file name


  /** <h3>A Single SHE Postprocessing Step</h3>

    The one- and two-dimensional meshes in other examples typically
    used several SHE iterations to obtain self-consistency of the
    Boltzmann equation and the Poisson equation. To keep
    the computational expense of this example under control,
    we assume that the electrostatic potential computed from the
    drift-diffusion model corresponds with the solution of the
    Boltzmann-Poisson system, hence we only compute a single
    solution of the Boltzmann equation using the SHE method
    for a given potential.

    Similar to the case of the drift-diffusion simulation above,
    we first need to set up the configuration.

    <h4>Prepare the SHE simulator configuration</h4>

    First we set up a new configuration object, enable
    electrons and holes, and specify that we want to use
    SHE for electons, but only a simple continuity equation for holes:
  **/
  std::cout << "* main(): Setting up first-order SHE (non-self consistent!)..." << std::endl;
  viennashe::config config;

  config.with_electrons(true); config.with_holes(true);
  config.set_electron_equation(viennashe::EQUATION_SHE); // Use SHE for electrons
  config.set_hole_equation(viennashe::EQUATION_CONTINUITY); // Use DD for holes

  /** We only want to run a first-order SHE simulation, hence we set the maximum expansion order accordingly: **/
  config.max_expansion_order(1);
  /** The energy spacing with respect to total energy H should be 31 meV,
      which amounts to half of the inelastic optical phonon scattering energy: **/
  config.energy_spacing(0.031 * viennashe::physics::constants::q);

  /** For the nonlinear solver we just set the maximum number of iteration to 1: **/
  config.nonlinear_solver().max_iters(1);

  /** <h4>Create and Run the SHE-Simulator</h4>
    The SHE simulator object is created in the same manner as the DD simulation object.
    The additional step here is to explicitly set the initial guesses:
    Quantities computed from the drift-diffusion simulation are passed to the SHE simulator object
    by means of the member function set_initial_guess().
    Then, the simulation is invoked using the member function run()
   **/
  std::cout << "* main(): Computing first-order SHE (requires about 5 GB RAM)..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);

  she_simulator.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());

  she_simulator.run();

  /** <h4>Write SHE Simulation Output</h4>

    Since this is a fully three-dimensional device simulation, we cannot
    write the full distribution function to VTK file(s), since VTK does not
    define four-dimensional cell types. Instead, one has to restrict
    to macroscopic values defined in the macroscopic domain.

    As for the drift-diffusion case, all macroscopic output quantities can be written
    to one or more VTK files using write_quantities_to_VTK_file():
  **/
  viennashe::io::write_quantities_to_VTK_file(she_simulator,            // simulator object
                                              "half-trigate_she_quan"); // file name

  /** If desired, one may also write individual quantities to VTK,
      for example the electrostatic potential or the electron density: **/
  viennashe::io::write_quantity_to_VTK_file(she_simulator.potential(),        device, "trigate_she_potential");
  viennashe::io::write_quantity_to_VTK_file(she_simulator.electron_density(), device, "trigate_she_electrons");

  /** Finally, print a small message to let the user know that everything succeeded **/
  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << std::endl;

  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;

}
