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
#include <stdexcept>

#include "tests/src/common.hpp"

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default mesh configuration:
#include "viennagrid/viennagrid.h"


/** \file simple_impurity_scattering.cpp Contains qulitative tests for impurity scattering.
 *  \test Qualitatively tests the impurity scattering by requiring certain physical relations to hold. Does not test absolute numbers!
 */

/** @brief Initalizes the device. Is typically modified by the user according to his/her needs.
*
* Can also be replaced by a reader that grabs all parameters from an external file
*
* @param device The device class that is to be initalized
* @param voltage The voltage to be applied on the right terminal; The other terminal is grounded
*/
template <typename DeviceType>
void init_device(DeviceType & device, double voltage)
{
  typedef typename DeviceType::mesh_type           MeshType;

  // STEP 2: Set doping
  std::cout << "* init_device(): Setting doping..." << std::endl;

  device.set_doping_n(1e28);
  device.set_doping_p(1e4);
  device.set_material(viennashe::materials::si());

  // STEP 3: Define contacts
  double gnd = 0.0;
  double vcc = voltage;// * 2.0;

  viennagrid_dimension cell_dim;
  viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

  viennagrid_element_id *cells_begin, *cells_end;
  viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

  if (cells_end - cells_begin < 3) throw std::runtime_error("The Mesh is too small. It contains less than 3 cells!");

  device.set_contact_potential(gnd, cells_begin[0]);
  device.set_material(viennashe::materials::metal(), cells_begin[0]);
  device.set_contact_potential(vcc, *(cells_end - 1));
  device.set_material(viennashe::materials::metal(), *(cells_end - 1));

}

int main()
{
  std::cout << viennashe::preamble() << std::endl;

  std::cout << "* main(): Initializing device..." << std::endl;
  viennashe::device device;

  //
  //  Generate device geometry using built-in device generator
  //
  {
    viennashe::util::device_generation_config generator_params;
    generator_params.add_segment(0.0,  1e-6, 51);  //start at x=0, length 1e-6, 51 points
    device.generate_mesh(generator_params);
  }

  //
  // Set up and initialize device
  //
  init_device(device, 0.1);

  //
  // Use a drift-diffusion simulation for obtaining an initial guess of the potential:
  //
  std::cout << "* main(): Creating DD simulator..." << std::endl;
  viennashe::config dd_cfg;
  dd_cfg.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);
  dd_cfg.nonlinear_solver().max_iters(50);
  dd_cfg.nonlinear_solver().damping(0.5);
  //dd_cfg.nonlinear_solver().set(viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);  //uncomment this to use Newton (might require stronger damping)

  viennashe::simulator dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;
  dd_simulator.run();

  //
  // A single SHE postprocessing step (fixed field)
  //
  std::cout << "* main(): Setting up SHE..." << std::endl;

  viennashe::config config;

  config.set_electron_equation(viennashe::EQUATION_SHE);
  config.with_electrons(true);

  config.with_holes(false);
  config.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  config.nonlinear_solver().max_iters(10);
  config.nonlinear_solver().damping(0.3);

  config.max_expansion_order(1);
  config.energy_spacing(0.0155 * viennashe::physics::constants::q);
  config.min_kinetic_energy_range(1.0 * viennashe::physics::constants::q);

  // configure scattering mechanisms. We use phonon scattering only for now.
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(false);
  config.scattering().impact_ionization().enabled(false);
  config.scattering().electron_electron(false);

  //config.dispersion_relation(viennashe::she::dispersion_relation_ids::ext_vecchi_dispersion);    //uncomment this to use the extended Vecchi model

  std::cout << "* main(): Computing SHE (electrons only) without ionized impurity and Modena model ..." << std::endl;
  viennashe::simulator she_simulator(device, config);

  she_simulator.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());

  she_simulator.run();

  std::cout << "* main(): Writing SHE result..." << std::endl;
  viennashe::io::gnuplot_edf_writer edfwriter;

  viennagrid_dimension cell_dim;
  viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

  viennagrid_element_id *cells_begin, *cells_end;
  viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

  double n_mid_without  = 0;
  double Jn_mid_without = 0;

  {
    viennashe::io::she_vtk_writer()(device,
                                    she_simulator.config(),
                                    she_simulator.quantities().electron_distribution_function(),
                                    "simple_impurity_scattering_edf1");

    edfwriter(device, viennashe::util::any_filter(), she_simulator.edf(viennashe::ELECTRON_TYPE_ID), "simple_impurity_scattering_edf1.dat");

    typedef viennashe::simulator::she_quantity_type she_quan_type;
    viennashe::she::current_density_wrapper<she_quan_type> Jn(device,
                  she_simulator.config(),
                  she_simulator.quantities().electron_distribution_function());

    viennashe::she::carrier_density_wrapper<she_quan_type> n(
                    she_simulator.config(),
                    she_simulator.quantities().electron_distribution_function());

    viennashe::io::write_cell_quantity_for_gnuplot(Jn, device, "simple_impurity_scattering_Jn1.dat");

    n_mid_without  = n(cells_begin[25])[0];
    Jn_mid_without = Jn(cells_begin[25])[0];
  }
  // ####################################################################################################################################

  config.scattering().ionized_impurity().enabled(true); // enable ionized impurity scattering

  std::cout << "* main(): Computing SHE (electrons only) WITH ionized impurity and Modena model ..." << std::endl;
  viennashe::simulator she_simulator_with(device, config);

  she_simulator_with.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator_with.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator_with.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());

  she_simulator_with.run();

  std::cout << "* main(): Writing SHE result..." << std::endl;

  viennashe::io::she_vtk_writer()(device,
                                  she_simulator_with.config(),
                                  she_simulator_with.quantities().electron_distribution_function(),
                                  "simple_impurity_scattering_edf2");

  //viennashe::io::gnuplot_edf_writer edfwriter;
  edfwriter(device, viennashe::util::any_filter(), she_simulator_with.edf(viennashe::ELECTRON_TYPE_ID), "simple_impurity_scattering_edf2.dat");

  typedef viennashe::simulator::she_quantity_type she_quan_type;
  viennashe::she::current_density_wrapper<she_quan_type> Jn(device,
                config,
                she_simulator_with.quantities().electron_distribution_function());

    viennashe::she::carrier_density_wrapper<she_quan_type> n(
                    config,
                    she_simulator_with.quantities().electron_distribution_function());

  viennashe::io::write_cell_quantity_for_gnuplot(Jn, device, "simple_impurity_scattering_Jn2.dat");

  double n_mid_with  = n(cells_begin[25])[0];
  double Jn_mid_with = Jn(cells_begin[25])[0];

  if(!viennashe::testing::fuzzy_equal(n_mid_without, n_mid_with))
  {
    std::cerr << "ERROR: Ionized impurity scattering should NOT change the electron density in a resistor!" << std::endl;
    std::cerr << "       The densities are not equal." << std::endl;
    std::cerr << n_mid_without << " != " << n_mid_with << std::endl;
    return EXIT_FAILURE;
  }

  if(viennashe::testing::fuzzy_equal(Jn_mid_without, Jn_mid_with, 1e-5))
  {
    std::cerr << "ERROR: Ionized impurity scattering should change the electron CURRENT density in a resistor!" << std::endl;
    std::cerr << "       It failed to do so." << std::endl;
    std::cerr << Jn_mid_without << " == " << Jn_mid_with << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cerr << "BUT: This is ok!" << std::endl;
  }

  if(Jn_mid_without < Jn_mid_with)
  {
    std::cerr << "ERROR: Ionized impurity scattering should DECREASE the electron CURRENT density in a resistor!" << std::endl;
    std::cerr << "       It failed to do so." << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}

