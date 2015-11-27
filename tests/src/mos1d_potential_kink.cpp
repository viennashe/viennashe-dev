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

// ViennaGrid default configuration and centroid() algorithm:
#include "viennagrid/viennagrid.h"


/** \file mos1d_potential_kink.cpp Contains tests to ensure that discont. permittivities are correctly accounted for in Poisson's equation
 *  \test Simulation of the potential in a purely 1d MOS to verify the kink due to discontinuous permittivities
 */


/** @brief A simple MOS 1D mesh generator for this particular test */
struct mos1d_mesh_generator
{

  mos1d_mesh_generator(double len_gate, double cs_gate, double len_oxide, double cs_ox, double len_bulk, double cs_bulk)
    : len_gate_(len_gate), cs_gate_(cs_gate), len_oxide_(len_oxide), cs_ox_(cs_ox), len_bulk_(len_bulk), cs_bulk_(cs_bulk)
  { }

  template<typename MeshT>
  void operator()(MeshT & mesh) const
  {
    viennashe::util::device_generation_config gconf;

    gconf.add_segment(0,                              len_gate_,  static_cast<unsigned long>(len_gate_ /cs_gate_ + 1.0) );
    gconf.add_segment(len_gate_,                      len_oxide_, static_cast<unsigned long>(len_oxide_/cs_ox_   + 1.0) );
    gconf.add_segment(len_gate_+len_oxide_,           len_bulk_,  static_cast<unsigned long>(len_bulk_ /cs_bulk_ + 1.0) );
    gconf.add_segment(len_gate_+len_oxide_+len_bulk_, 1e-9, 10 );

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

/** @brief Initalizes the MOS 1D device for testing. The Bulk-Contact is grounded.
*
* @param device The device class that is to be initalized
* @param Vg_init The initial gate potential w.r.t. ground (bulk-contact)
*/
template <typename DeviceType>
void init_device(DeviceType & device, double Vg_init)
{
  typedef typename DeviceType::segment_type        SegmentType;

  SegmentType const & gate    = device.segment(0);
  SegmentType const & oxide   = device.segment(1);
  SegmentType const & silicon = device.segment(2);
  SegmentType const & bulk    = device.segment(3);

  std::cout << "* init_device(): Setting material ..." << std::endl;

  device.set_material(viennashe::materials::metal(), gate);
  device.set_material(viennashe::materials::metal(), bulk);
  device.set_material(viennashe::materials::hfo2(),  oxide);
  device.set_material(viennashe::materials::si(),    silicon);

  std::cout << "* init_device(): Setting doping (per cell) ..." << std::endl;

  //intrinsic doping for silicon to compensate built-in potentials
  device.set_doping_n(1e16, silicon);
  device.set_doping_p(1e16, silicon);

  std::cout << "* init_device(): Setting contact potentials (per cell) ..." << std::endl;

  device.set_contact_potential(0.0,      bulk);
  device.set_contact_potential(Vg_init,  gate);

  std::cout << "* init_device(): DONE!" << std::endl;

} // init_device()

int main()
{
  typedef viennashe::device<viennagrid_mesh>     DeviceType;

  //typedef viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;
  //typedef viennagrid::result_of::iterator<CellContainer>::type        CellIterator;
  //typedef viennagrid::result_of::const_facet_range<MeshType>::type    FacetContainer;

  std::cout << "* main(): Initializing device ..." << std::endl;

  DeviceType device;

  const double len_gate  = 1e-9;
  const double cs_gate   = 1e-10;
  const double len_oxide = 1e-9;
  const double cs_ox     = 2e-10;
  const double len_bulk  = 1e-9;
  const double cs_bulk   = 1e-10;

  double Vg = 0.3;

  mos1d_mesh_generator mosgen(len_gate, cs_gate, len_oxide, cs_ox, len_bulk, cs_bulk);
  device.generate_mesh(mosgen);

  init_device(device, Vg);

  //
  // Use a drift-diffusion simulation for obtaining an initial guess of the potential:
  //
  std::cout << "* main(): Creating DD simulator..." << std::endl;

  viennashe::config dd_cfg;
  dd_cfg.with_electrons(false);
  dd_cfg.with_holes(false);
  dd_cfg.nonlinear_solver().max_iters(20);
  dd_cfg.nonlinear_solver().damping(1.0);

  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;

  dd_simulator.run();

  viennagrid_dimension cell_dim;
  viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

  viennagrid_element_id *cells_begin, *cells_end;
  viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

  double centroid13[3], centroid12[3];
  viennagrid_element_centroid(device.mesh(), cells_begin[12], centroid12);
  viennagrid_element_centroid(device.mesh(), cells_begin[13], centroid13);
  double pot_gate_left  = dd_simulator.potential().at(cells_begin[12]);
  double pot_gate_right = dd_simulator.potential().at(cells_begin[13]);
  double dx_gate = std::fabs(centroid13[0] - centroid12[0]);
  double field_gate = (pot_gate_left - pot_gate_right) / dx_gate;

  double centroid23[3], centroid22[3];
  viennagrid_element_centroid(device.mesh(), cells_begin[22], centroid22);
  viennagrid_element_centroid(device.mesh(), cells_begin[23], centroid23);
  double pot_bulk_left  = dd_simulator.potential().at(cells_begin[22]);
  double pot_bulk_right = dd_simulator.potential().at(cells_begin[23]);
  double dx_bulk = std::fabs(centroid23[0] - centroid22[0]);
  double field_bulk = (pot_bulk_left - pot_bulk_right) / dx_bulk;


  std::cout << "Field in gate: " << field_gate << std::endl;
  std::cout << "Field in bulk: " << field_bulk << std::endl;
  std::cout << "Ratio of electric fields: " << field_gate / field_bulk << std::endl;
  double permittivity_ox = viennashe::materials::permittivity(device.get_material(cells_begin[13]));
  double permittivity_bulk = viennashe::materials::permittivity(device.get_material(cells_begin[23]));
  std::cout << "Ratio of permittivities: " << permittivity_bulk / permittivity_ox  << std::endl;

  viennashe::io::gnuplot_writer gpwriter;
  gpwriter(device, viennashe::util::any_filter(), dd_simulator.potential(), "mos1d_potential_kink.dat");

  std::cout << std::endl;
  if ( std::abs(field_gate / field_bulk - permittivity_bulk / permittivity_ox) > 0.01 )
  {
    std::cout << "TEST FAILED: Ratios of fields do not correspond to ratio of permittivities!" << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "TEST PASSED: Ratios match!" << std::endl;
  std::cout << std::endl;

  std::cout << "*********************************************************" << std::endl;
  std::cout << "*               Test finished successfully              *" << std::endl;
  std::cout << "*********************************************************" << std::endl;

  return EXIT_SUCCESS;
}
