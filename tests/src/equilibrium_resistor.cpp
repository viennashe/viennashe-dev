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

#include "viennashe/forwards.h"

#include "tests/src/common.hpp"

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default configurations and centroid() algorithm:
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/algorithm/centroid.hpp"


/** \file equilibrium_resistor.cpp Contains a test of the SHE of the BTE in a resistor at equilibrium
 *  \test Test of the SHE of the BTE in a resistor at equilibrium. The EDF must be Maxwellian
 */

/**
 * @brief A simple wrapper for the kinetic energy of charge carriers
 * @param q The SHE quantity having a member called 'kinetic_energy(element, index_H, carrier_type)
 * @param ctype The carrier type for which to wrap the kinetic energy
 */
template <typename SHEQuanT>
class kinetic_energy_wrapper
{
  public:
    kinetic_energy_wrapper(SHEQuanT const & q, viennashe::carrier_type_id ctype) : quan_(q), carrier_type_(ctype) {}

    template <typename ElementType>
    double operator()(ElementType const & el, std::size_t index_H) const
    {
      return quan_.kinetic_energy(el, index_H, carrier_type_);
    }

  private:
    SHEQuanT const & quan_;
    viennashe::carrier_type_id carrier_type_;
};


/** @brief Initalizes the device. Is typically modified by the user according to his/her needs.
*
* Can also be replaced by a reader that grabs all parameters from an external file
*
* @param device The device class that is to be initalized
* @param len_x  The length of the resistor in meter
*/
template <typename DeviceType>
void init_device(DeviceType & device, double len_x)
{
  typedef typename DeviceType::mesh_type           MeshType;

  // STEP 2: Set doping
  std::cout << "* init_device(): Setting doping..." << std::endl;

  // STEP 1: Set material properties:
  device.set_doping_n(1e16);
  device.set_doping_p(1e16);
  device.set_material(viennashe::materials::si());


  // STEP 2: Define contacts
  typedef typename viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

  double gnd = 0.0;
  double vcc = 0.0;

  CellContainer cells(device.mesh());
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    //left contact:
    if (viennagrid::centroid(*cit)[0] < 0.1 * len_x)
    {
      device.set_contact_potential(gnd, *cit);
    }

    //right contact:
    if (viennagrid::centroid(*cit)[0] > 0.9 * len_x)
      device.set_contact_potential(vcc, *cit);
  }

}

/**
 * @brief Tests whether the EDF at a certain given point is Maxwellian
 * @param carrier_type The carrier type for which to test the EDF
 * @param device The device
 * @param quan The SHE quantity, that is the EDF
 * @param conf The configuration used for simulation
 * @param cell The cell for which to test (we use the centroid of the cell)
 * @param index_H The H index at which to perform the test
 * @return True if the EDF is okay, that is maxwellian, else false
 */
template <typename DeviceType, typename SHEQuantity, typename CellType>
int test_result_at_point(viennashe::carrier_type_id carrier_type,
                         DeviceType const & device,
                         SHEQuantity const & quan,
                         viennashe::config const & conf,
                         CellType const & cell,
                         std::size_t index_H)
{
  const double T   = device.get_lattice_temperature(cell);
  const double tol = 1e-3;
  double kinetic_energy = quan.get_kinetic_energy(cell, index_H);
  viennashe::math::SphericalHarmonic Y_00(0,0);

  double ref = std::exp(- kinetic_energy / viennashe::physics::constants::kB / T)
                * quan.get_values(cell, 1)[0]; //for normalization

  // Test 1: Check Maxwell distribution of zeroth-order SHE coefficients:
  if ( !viennashe::testing::fuzzy_equal(quan.get_values(cell, index_H)[0], ref, tol)
       || (quan.get_values(cell, index_H)[0] < quan.get_values(cell, index_H)[0])
       || (quan.get_values(cell, index_H)[0] > quan.get_values(cell, index_H)[0])
       )
  {
    std::cerr << "* ERROR: Maxwell distribution test failed at vertex " << cell << " with index_H = " << index_H << std::endl;
    return EXIT_FAILURE;
  }

  // Test 2: Interpolated first-order expansion coefficients should be to zero (up to round-off):
  viennashe::she::interpolated_she_df_wrapper<DeviceType, SHEQuantity> interpolated_she_df(device, conf, quan);
  for (long m = -1; m <= 1; ++m)
  {
    double f_1m = interpolated_she_df(cell, kinetic_energy, 1, m, index_H);
    if ( std::abs(f_1m) > 1e-18 )
    {
      std::cerr << "* ERROR: Check for vanishing first-order failed at vertex " << cell << " with index_H = " << index_H << std::endl;
      std::cerr << " value: " << f_1m << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Test 3: Generalized distribution function evaluated at random angles should be the same (up to round-off)
  viennashe::she::generalized_df_wrapper<DeviceType, SHEQuantity> generalized_df(device, conf, quan);
  double gdf_val1 = generalized_df(cell, kinetic_energy, 0.123, 1.02);
  double gdf_val2 = generalized_df(cell, kinetic_energy, 0.420, 3.60);
  double gdf_val3 = generalized_df(cell, kinetic_energy, 0.360, 4.02);
  ref *= conf.dispersion_relation(carrier_type).density_of_states(kinetic_energy) * Y_00(0, 0);
  if (   ! viennashe::testing::fuzzy_equal(gdf_val1, ref, tol)
      || ! viennashe::testing::fuzzy_equal(gdf_val2, ref, tol)
      || ! viennashe::testing::fuzzy_equal(gdf_val3, ref, tol)
     )
  {
    std::cerr << "* ERROR: Check for generalized distribution function failed at vertex " << cell << " with index_H = " << index_H << std::endl;
    std::cerr << " val1: " << gdf_val1 << std::endl;
    std::cerr << " val2: " << gdf_val2 << std::endl;
    std::cerr << " val3: " << gdf_val3 << std::endl;
    std::cerr << "  ref: " <<      ref << std::endl;
    return EXIT_FAILURE;
  }

  // Test 4: Now check the same for the generalized energy distribution function wrapper:
  viennashe::she::generalized_edf_wrapper<DeviceType, SHEQuantity> generalized_edf(conf, quan);
  double edf_val1 = generalized_edf(cell, kinetic_energy);
  double edf_val2 = generalized_edf(cell, kinetic_energy);
  double edf_val3 = generalized_edf(cell, kinetic_energy);
  if (   ! viennashe::testing::fuzzy_equal(edf_val1, ref, tol)
      || ! viennashe::testing::fuzzy_equal(edf_val2, ref, tol)
      || ! viennashe::testing::fuzzy_equal(edf_val3, ref, tol)
     )
  {
    std::cerr << "* ERROR: Check for generalized energy distribution function failed at vertex " << cell << " with index_H = " << index_H << std::endl;
    std::cerr << " val1: " << edf_val1 << std::endl;
    std::cerr << " val2: " << edf_val2 << std::endl;
    std::cerr << " val3: " << edf_val3 << std::endl;
    std::cerr << "  ref: " <<      ref << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}


inline int simulate(double temperature)
{
  typedef viennagrid::quadrilateral_2d_mesh                     MeshType;
  typedef viennashe::device<MeshType>                           DeviceType;

  std::cout << "* main(): Creating device..." << std::endl;
  DeviceType device;

  // STEP 1: Generate device:
  viennashe::util::device_generation_config generator_params;
  generator_params.add_segment(0.0, 6.0e-7, 10,   //start at x=, length, points
                               0.0, 6.0e-7,  2);  //start at y=, length, points
  device.generate_mesh(generator_params);
  device.set_lattice_temperature(temperature);

  std::cout << "* main(): Initializing device..." << std::endl;
  init_device(device, generator_params.at(0).get_length_x());

  std::cout << "* main(): Creating DD simulator..." << std::endl;
  viennashe::config dd_cfg;
  dd_cfg.with_electrons(true);
  dd_cfg.set_electron_equation(viennashe::EQUATION_CONTINUITY);

  dd_cfg.with_holes(true);
  dd_cfg.set_hole_equation(viennashe::EQUATION_CONTINUITY);
  //dd_cfg.gummel_iters(50);
  viennashe::simulator<DeviceType> dd_simulator(device, dd_cfg);

  std::cout << "* main(): Launching simulator..." << std::endl;
  dd_simulator.run();

  viennashe::io::write_quantities_to_VTK_file(dd_simulator, "equilibrium-resistor_dd_quan");

  // A single SHE solution step:
  std::cout << "* main(): Setting up SHE..." << std::endl;
  //size_t energy_points = 10;
  viennashe::config config;
  //config.max_expansion_order(3);
  config.with_electrons(true);
  config.with_holes(true);
  config.set_electron_equation(viennashe::EQUATION_SHE);
  config.set_hole_equation(viennashe::EQUATION_SHE);
  config.scattering().ionized_impurity().enabled(false);
  //config.energy_levels(energy_points);
  config.energy_spacing(6.2 * viennashe::physics::constants::q / 1000.0);
  config.nonlinear_solver().max_iters(1);

  std::cout << "* main(): Computing SHE..." << std::endl;
  viennashe::simulator<DeviceType> she_simulator(device, config);
  she_simulator.set_initial_guess(viennashe::quantity::potential(),        dd_simulator.potential());
  she_simulator.set_initial_guess(viennashe::quantity::electron_density(), dd_simulator.electron_density());
  she_simulator.set_initial_guess(viennashe::quantity::hole_density(),     dd_simulator.hole_density());
  she_simulator.run();

  std::cout << "Energy range: [" << she_simulator.quantities().electron_distribution_function().get_value_H(0) << ", "
                                 << she_simulator.quantities().electron_distribution_function().get_value_H(she_simulator.quantities().electron_distribution_function().get_value_H_size() - 1)
                                 << "] (" << she_simulator.quantities().electron_distribution_function().get_value_H_size() << " discrete energies)" << std::endl;

  std::cout << "* main(): Writing SHE result..." << std::endl;
  std::stringstream ss;
  ss << "equilibrium_resistor_" << temperature << "_edf";
  viennashe::io::she_vtk_writer<DeviceType>()(device,
                                              she_simulator.config(),
                                              she_simulator.quantities().electron_distribution_function(),
                                              ss.str());

  //
  // check result: Must equal kinetic energy:
  //
  typedef viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

  //VectorType she_result = she_simulator.edf().vector();
  //std::cout << "she_result: " << she_result << std::endl;

  CellContainer cells(device.mesh());
  for (CellIterator cit = cells.begin();
       cit != cells.end();
       ++cit)
  {
    for (size_t index_H=0; index_H < she_simulator.quantities().electron_distribution_function().get_value_H_size(); ++index_H)
    {
      long index = she_simulator.quantities().electron_distribution_function().get_unknown_index(*cit, index_H);// indices_on_vertex[index_H];
      if (index > -1)
      {
        if (she_simulator.config().with_electrons())
        {
          if (test_result_at_point(viennashe::ELECTRON_TYPE_ID, device, she_simulator.quantities().electron_distribution_function(), config, *cit, index_H) != EXIT_SUCCESS)
            return EXIT_FAILURE;
        }
        else if (she_simulator.config().with_holes())
        {
          if (test_result_at_point(viennashe::HOLE_TYPE_ID, device, she_simulator.quantities().hole_distribution_function(), config, *cit, index_H) != EXIT_SUCCESS)
            return EXIT_FAILURE;
        }
      }
    }
    std::cout << "Tests passed for cell " << cit->id().get() << std::endl;
  }

  return EXIT_SUCCESS;
}


int main()
{
  // Simulate at 300 K and at 400 K to make sure temperature is considered accordingly:
  std::cout << "--------- Testing at 300 K --------------" << std::endl;
  int retval = simulate(300.0);

  if (retval != EXIT_SUCCESS)
    return retval;

  std::cout << std::endl;
  std::cout << "--------- Testing at 400 K --------------" << std::endl;
  retval = simulate(400.0);
  if (retval != EXIT_SUCCESS)
    return retval;

  std::cout << "*******************************" << std::endl;
  std::cout << "* Test finished successfully! *" << std::endl;
  std::cout << "*******************************" << std::endl;

  return EXIT_SUCCESS;
}

