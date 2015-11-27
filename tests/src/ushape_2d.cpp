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

// ViennaGrid default configurations:
#include "viennagrid/viennagrid.h"

/** \file ushape_2d.cpp Tests charge conservation in a U-shaped configuration.
 *  \test Tests charge conservation on a block of silicon with two contacts. The current is expected to flow in a U-shaped config
 */



/** @brief Initalizes the device. Is typically modified by the user according to his/her needs.
*
* Can also be replaced by a reader that grabs all parameters from an external file
*
* @param device The device class that is to be initalized
*/
template <typename DeviceType>
void init_device(DeviceType & device)
{
  typedef typename DeviceType::segment_type          SegmentType;

  SegmentType contact_left  = device.segment(1);
  SegmentType oxide         = device.segment(2);
  SegmentType contact_right = device.segment(3);
  SegmentType body          = device.segment(4);

  device.set_material(viennashe::materials::si(), body);

  device.set_material(viennashe::materials::hfo2(), oxide);

  device.set_material(viennashe::materials::metal(), contact_left);
  device.set_material(viennashe::materials::metal(), contact_right);

  device.set_doping_n(1e24, body);
  device.set_doping_p(1e8,  body);


  // Set contact potentials
  device.set_contact_potential(0.0, contact_left);
  device.set_contact_potential(0.5, contact_right);

}

int main()
{
  typedef viennashe::device<viennagrid_mesh>           DeviceType;
  typedef DeviceType::segment_type                     SegmentType;

  std::cout << viennashe::preamble() << std::endl;

  std::cout << "* main(): Creating device..." << std::endl;
  DeviceType device;
  device.load_mesh("../../tests/data/ushape2d/ushape125.mesh");

  viennagrid_element_id *vertices_begin, *vertices_end;
  viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

  viennagrid_element_id *edges_begin, *edges_end;
  viennagrid_mesh_elements_get(device.mesh(), 1, &edges_begin, &edges_end);

  std::cout << "Vertices: " << vertices_end - vertices_begin << std::endl;
  std::cout << "Edges:    " << edges_end - edges_begin       << std::endl;

  // scale device from meter to nanometer
  device.scale(1e-6);

  init_device(device);

  //
  // Use a drift-diffusion simulation for obtaining an initial guess of the potential:
  //
  std::cout << "* main(): Creating DD simulator..." << std::endl;
  viennashe::config dd_cfg;
  dd_cfg.with_electrons(true);
  dd_cfg.with_holes(true);
  dd_cfg.nonlinear_solver().max_iters(300);
  dd_cfg.nonlinear_solver().damping(0.4);
  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching DD simulator..." << std::endl;
  dd_simulator.run();

  //
  // Write drift-diffusion quantities to file for inspection (e.g. using ParaView)
  //
  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "ushape2d_dd_quan");

  // Check current conservation:
  std::cout << "#" << std::endl;
  std::cout << "# Checking electron current: " << std::endl;
  std::cout << "#" << std::endl;
  viennashe::check_current_conservation(device,
                                        viennashe::ELECTRON_TYPE_ID,
                                        dd_simulator.potential(),
                                        dd_simulator.electron_density(),
                                        viennashe::models::create_constant_mobility_model(device, 0.1430));

  std::cout << "#" << std::endl;
  std::cout << "# Checking hole current: " << std::endl;
  std::cout << "#" << std::endl;
  viennashe::check_current_conservation(device,
                                        viennashe::HOLE_TYPE_ID,
                                        dd_simulator.potential(),
                                        dd_simulator.hole_density(),
                                        viennashe::models::create_constant_mobility_model(device, 0.0460));


  //
  // Self-consistent SHE simulations
  //

  std::cout << "* main(): Setting up SHE..." << std::endl;

  //std::size_t energy_points = 100;
  viennashe::config config;
  config.set_electron_equation(viennashe::EQUATION_SHE);
  config.with_electrons(true);
  config.set_hole_equation(viennashe::EQUATION_CONTINUITY);
  config.with_holes(true);
  config.nonlinear_solver().max_iters(40);   //Use higher iteration counts as needed
  config.nonlinear_solver().damping(0.4);
  config.max_expansion_order(1);
  //config.she_boundary_conditions().type(viennashe::BOUNDARY_DIRICHLET);
  //config.energy_levels(energy_points);
  config.energy_spacing(31.0 * viennashe::physics::constants::q / 1000.0);  //energy spacing of 31 meV
  config.linear_solver().set(viennashe::solvers::linear_solver_ids::serial_linear_solver);        //a serial linear solver
  //config.linear_solver().set(viennashe::solvers::linear_solver_ids::parallel_linear_solver);     //use this if OpenMP is available
  //config.linear_solver().set(viennashe::solvers::linear_solver_ids::gpu_parallel_linear_solver); //use this for GPU support

  //
  // Instantiate and launch the SHE simulator
  //
  std::cout << "* main(): Computing SHE..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);
  she_simulator.set_initial_guess(viennashe::quantity::potential(), dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(), dd_simulator.hole_density());
  she_simulator.run();

  //
  // Write output
  //
  std::cout << "* main(): Writing SHE result..." << std::endl;

  // Writes EDF to VTK file
  viennashe::io::she_vtk_writer<DeviceType>()(device,
                                              she_simulator.config(),
                                              she_simulator.quantities().electron_distribution_function(),
                                              "ushape2d_edf");

  viennashe::io::write_quantities_to_VTK_file(she_simulator, "ushape2d_she_quan");


  SegmentType const & contact_left  = device.segment(1);
  SegmentType const & contact_right = device.segment(3);
  SegmentType const & semiconductor = device.segment(4);

  double current_left = viennashe::get_terminal_current(device,
                                                        she_simulator.config(),
                                                        she_simulator.quantities().electron_distribution_function(),
                                                        semiconductor,
                                                        contact_left);
  double current_right = viennashe::get_terminal_current(device,
                                                         she_simulator.config(),
                                                         she_simulator.quantities().electron_distribution_function(),
                                                         semiconductor,
                                                         contact_right);

  std::cout << "Current left:  " << current_left  << std::endl;
  std::cout << "Current right: " << current_right << std::endl;

  viennashe::she::check_current_conservation(device,
                                             she_simulator.config(),
                                             she_simulator.quantities().electron_distribution_function());

  if ( std::fabs(current_left + current_right) / std::max(std::fabs(current_left), std::fabs(current_right)) > 1e-5)
  {
    std::cerr << "CURRENT MISMATCH DETECTED" << std::endl;
    std::cerr << "Norm of error: " << std::fabs(current_left - current_right) / std::max(std::fabs(current_left), std::fabs(current_right)) << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "SUCCESS: CURRENT MATCHES!" << std::endl;

  std::cout << "* main(): Results can now be viewed with your favorite VTK viewer (e.g. ParaView)." << std::endl;
  std::cout << "* main(): Don't forget to scale the z-axis by about a factor of 1e11 when examining the distribution function." << std::endl;
  std::cout << std::endl;
  std::cout << "****************************************************" << std::endl;
  std::cout << "*           Test finished successfully             *" << std::endl;
  std::cout << "****************************************************" << std::endl;

  return EXIT_SUCCESS;
}

