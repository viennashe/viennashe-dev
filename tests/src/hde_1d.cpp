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
#include <cmath>
#include <vector>

#include "tests/src/common.hpp"
#include "mos1d_dg.hpp"

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default configurations:
#include "viennagrid/viennagrid.h"


/** \file hde_1d.cpp Contains a simple test of the heat diffusion equation (essentially a poisson equation).
 *  \test A very simple test of the HDE implementation, just to make sure it does not die.
 */


/**
 * @brief A Simple test mesh generator... generates a 1D MOS structure for the HDE test
 */
struct hde_test_mesh_generator
{

  hde_test_mesh_generator(double len_gate, double cs_gate, double len_oxide, double cs_ox, double len_bulk, double cs_bulk)
    : len_gate_(len_gate), cs_gate_(cs_gate), len_oxide_(len_oxide), cs_ox_(cs_ox), len_bulk_(len_bulk), cs_bulk_(cs_bulk)
  { }

  template<typename MeshT>
  void operator()(MeshT & mesh) const
  {
    viennashe::util::device_generation_config gconf;

    gconf.add_segment(0,                              len_gate_,        static_cast<unsigned long>(std::ceil(len_gate_/cs_gate_)+1) );
    gconf.add_segment(len_gate_,                      len_oxide_,       static_cast<unsigned long>(std::ceil(len_oxide_/cs_ox_))    );
    gconf.add_segment(len_gate_+len_oxide_,           len_bulk_,        static_cast<unsigned long>(std::ceil(len_bulk_/cs_bulk_))   );
    gconf.add_segment(len_gate_+len_oxide_+len_bulk_, len_bulk_ + 2e-9, static_cast<unsigned long>(std::ceil(2e-9/cs_bulk_)+1)      );

    viennashe::util::generate_device(mesh, gconf);
  }

private:
  double len_gate_;
  double cs_gate_;
  double len_oxide_;
  double cs_ox_;
  double len_bulk_;
  double cs_bulk_;

};

/**
 * @brief Initializes the device for a HDE test, where the left and right
 *        terminals are supposed to be at different lattice temperatures
 * @param device A reference to the device
 * @param TL_left The temperature of the left terminal
 * @param TL_right The temperature of the right terminal
 */
template <typename DeviceType>
void init_device(DeviceType & device, double TL_left, double TL_right)
{
  typedef typename DeviceType::mesh_type           MeshType;
  typedef typename DeviceType::segment_type        SegmentType;

  SegmentType gate    = device.segment(0);
  SegmentType oxide   = device.segment(1);
  SegmentType silicon = device.segment(2);
  SegmentType bulk    = device.segment(3);

  std::cout << "* init_device(): Setting material ..." << std::endl;

  device.set_material(viennashe::materials::metal(), gate);
  device.set_material(viennashe::materials::metal(), bulk);
  device.set_material(viennashe::materials::sio2(),  oxide);
  device.set_material(viennashe::materials::si(),    silicon);

  std::cout << "* init_device(): Setting doping (per cell) ..." << std::endl;

  device.set_doping_n(1e8,  silicon);
  device.set_doping_p(1e23, silicon);

  std::cout << "* init_device(): Setting contact potentials (per cell) ..." << std::endl;

  device.set_contact_potential(0.0,  bulk);
  device.set_contact_potential(0.0,  gate);

  std::cout << "* init_device(): Setting temperature boundary conditions (per cell) ..." << std::endl;

  viennagrid_dimension cell_dim;
  viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

  viennagrid_element_id *cells_begin, *cells_end;
  viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);
  for (viennagrid_element_id *cit  = cells_begin;
                              cit != cells_end;
                            ++cit)
  {
    device.set_lattice_temperature(300.0, *cit);
  }

  device.set_lattice_temperature(TL_left,  gate);
  device.set_lattice_temperature(TL_right, bulk);

  std::cout << "* init_device(): DONE!" << std::endl;

} // init_device()

int main()
{
  typedef viennashe::device<viennagrid_mesh>    DeviceType;

  std::cout << viennashe::preamble() << std::endl;

  std::cout << "* main(): Creating mesh ..." << std::endl;

  DeviceType device;

  const double len_gate  = 1e-9;
  const double cs_gate   = 0.01e-9;
  const double len_oxide = 1e-9;
  const double cs_ox     = 0.01e-9;
  const double len_bulk  = 100e-9;
  const double cs_bulk   = 2e-9;

  hde_test_mesh_generator mosgen(len_gate, cs_gate, len_oxide, cs_ox, len_bulk, cs_bulk);
  device.generate_mesh(mosgen);

  std::cout << "* main(): Initializing device ..." << std::endl;

  init_device(device, 350.0, 300.0);

  std::cout << "* main(): Simulation setup ..." << std::endl;

  viennashe::config dd_cfg;
  dd_cfg.with_holes(true);
  dd_cfg.with_electrons(true);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);
  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);

  dd_cfg.with_hde(true); // USE THE HDE !

  // We need the dense solver here ...
  dd_cfg.linear_solver().set(viennashe::solvers::linear_solver_ids::dense_linear_solver);
  dd_cfg.nonlinear_solver().max_iters(30);
  dd_cfg.nonlinear_solver().damping(0.8);

  viennashe::simulator<DeviceType> simulator(device, dd_cfg);

  std::cout << "* main(): Running simulation ..." << std::endl;

  simulator.run();

  viennashe::io::write_quantities_to_VTK_file(simulator, "hde_1d_test_quan");
  viennashe::io::write_cell_quantity_for_gnuplot(simulator.quantities().lattice_temperature(),  device, "hde_1d_temperature.dat" );

  std::cout << "* main(): Testing ..." << std::endl;

  bool ok = test_result(simulator.quantities().lattice_temperature(),  device, "../../tests/data/hde/simple_mos.dat"  );
  if(!ok) return EXIT_FAILURE;

  std::cout << "* main(): Tests OK!" << std::endl;


  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*           ViennaSHE finished successfully             *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}

