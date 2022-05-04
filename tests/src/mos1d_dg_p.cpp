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

#include "tests/src/mos1d_dg.hpp"

/** \file mos1d_dg_p.cpp Contains tests of a purely 1d MOS (metal-SiO2-Si-metal) using the density gradient model for quantum corrections.
 *  \test This is a test for the density gradient model. The test is performed for holes in 1D MOS with DD only.
 */

int main()
{
  typedef viennagrid::line_1d_mesh                               MeshType;

  typedef viennashe::device<MeshType>     DeviceType;

  typedef viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type        CellIterator;
  typedef viennagrid::result_of::const_facet_range<MeshType>::type    FacetContainer;

  std::cout << viennashe::preamble() << std::endl;

  std::cout << "* main(): Creating mesh ..." << std::endl;

  DeviceType device;

  const double len_gate  = 1e-9;
  const double cs_gate   = 0.01e-9;
  const double len_oxide = 1e-9;
  const double cs_ox     = 0.05e-9;
  const double len_bulk  = 100e-9;
  const double cs_bulk   = 1e-9;

  double Vg = -0.2;

  bool ok   = false;

  mos1d_mesh_generator mosgen(len_gate, cs_gate, len_oxide, cs_ox, len_bulk, cs_bulk);
  device.generate_mesh(mosgen);

  //                       ND     NA
  init_device(device, Vg, 3e23,  1e8 );

  //
  // DEBUG OUTPUT
  std::cout << std::endl;
  CellContainer cells(device.mesh());
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    std::cout << cit->id() << " => " << device.get_material(*cit)
              << " ## Nd = " << device.get_doping_n(*cit)
              << " ## Na = " << device.get_doping_p(*cit)
              << " ## has bnd cond: " << (device.has_contact_potential(*cit)== true ? "TRUE" : "FALSE")
              << std::endl;
  } // for cells
  FacetContainer facets(device.mesh());
  for (std::size_t i = 0; i < facets.size(); ++i)
    std::cout << facets[i] << std::endl;

  //
  // Use a drift-diffusion simulation
  std::cout << "* main(): Creating DD simulator..." << std::endl;

  viennashe::config dd_cfg;
  dd_cfg.with_holes(true);
  dd_cfg.with_electrons(true);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  // We need the dense solver here ...
  dd_cfg.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);
  dd_cfg.nonlinear_solver().max_iters(30);
  dd_cfg.nonlinear_solver().damping(0.6);

  dd_cfg.quantum_correction(true);      // solve density gradient equations
  dd_cfg.with_quantum_correction(true); // apply correction potentials from density gradient

  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;

  dd_simulator.run();

  viennashe::io::gnuplot_writer gpwriter;
  gpwriter(device, viennashe::util::any_filter(), dd_simulator.electron_density(), "mos1d_dg_n_dd_p.dat");
  gpwriter(device, viennashe::util::any_filter(), dd_simulator.hole_density(),     "mos1d_dg_n_dd_p.dat");

  std::cout << "* main(): Testing results (DD) ..." << std::endl;

  ok = test_result(dd_simulator.dg_pot_p(),  device, "../../tests/data/dg/p/dd_dgpot_holes.dat"     );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.dg_pot_n(),  device, "../../tests/data/dg/p/dd_dgpot_electrons.dat" );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.potential(), device, "../../tests/data/dg/p/mos1d_dd_pot.dat"       );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.electron_density(), device, "../../tests/data/dg/p/mos1d_dd_electrons.dat");
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.hole_density(),     device, "../../tests/data/dg/p/mos1d_dd_holes.dat");
  if(!ok) return EXIT_FAILURE;


#if 0 // DISABLED DUE TO CONVERGENCY PROBLEMS

  //
  // Start self-consistent SHE simulations
  //
  std::cout << "* main(): Setting up SHE..." << std::endl;

  viennashe::config config;
  config.with_electrons(true);
  config.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  config.with_holes(true);
  config.set_hole_equation(viennashe::EQUATION_SHE);
  config.with_traps(false);
  config.scattering().acoustic_phonon().enabled(true);
  config.scattering().optical_phonon().enabled(true);
  config.scattering().ionized_impurity().enabled(true);
  config.scattering().impact_ionization().enabled(false);
  config.scattering().electron_electron(false);

  //config.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);
  config.nonlinear_solver().max_iters(30);
  config.nonlinear_solver().damping(0.5);     //use smaller damping with Newton
  config.max_expansion_order(1);
  config.linear_solver().max_iters(2000);
  config.linear_solver().ilut_drop_tolerance(1e-2);
  config.min_kinetic_energy_range(viennashe::physics::constants::q * 0.5);
  config.energy_spacing(viennashe::physics::constants::q * 0.031);  //energy spacing of 31meV

  config.quantum_correction(true);      // calculate a correction potential
  config.with_quantum_correction(true); // use the correction potential for the transport model

  //
  // Launch solver
  //
  std::cout << "* main(): Computing SHE..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);
  she_simulator.set_initial_guess(viennashe::quantity::potential(), dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(), dd_simulator.hole_density());
  //she_simulator.set_initial_guess(viennashe::quantity::density_gradient_electron_correction(), dd_simulator.dg_pot_n());
  //she_simulator.set_initial_guess(viennashe::quantity::density_gradient_hole_correction(), dd_simulator.dg_pot_p());
  she_simulator.run();

  //
  // TEST
  //
  std::cout << "* main(): Testing results (SHE) ..." << std::endl;

  ok = test_result(dd_simulator.dg_pot_p(),  device, "../../tests/data/dg/p/she_dgpot_holes.dat"     );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.dg_pot_n(),  device, "../../tests/data/dg/p/she_dgpot_electrons.dat" );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.potential(), device, "../../tests/data/dg/p/mos1d_she_pot.dat"       );
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.electron_density(), device, "../../tests/data/dg/p/mos1d_she_electrons.dat");
  if(!ok) return EXIT_FAILURE;
  ok = test_result(dd_simulator.hole_density(),     device, "../../tests/data/dg/p/mos1d_she_holes.dat");
  if(!ok) return EXIT_FAILURE;
#endif

  std::cout << "... \\o/ SUCCESS \\o/ ..." << std::endl;

  return EXIT_SUCCESS;
}
