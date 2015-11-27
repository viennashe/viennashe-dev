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

// C++ includes
#include "viennashe_all.hpp"

// C includes
#include "libviennashe/include/simulator.h"


namespace libviennashe
{

  /**
   * @brief Set the solution of one simulator as initial guess for another simulator
   * @param sim   The simulator, which gets its initial guesses set
   * @param other_sim The simulator, providing the solution
   */
  template < typename SimulatorT >
  void set_initial_guess(SimulatorT & sim, SimulatorT const & other_sim)
  {
    sim.set_initial_guess(viennashe::quantity::potential(),        other_sim.potential());
    sim.set_initial_guess(viennashe::quantity::electron_density(), other_sim.electron_density());
    sim.set_initial_guess(viennashe::quantity::hole_density(),     other_sim.hole_density());
  }

  /**
   * @brief Sets the inital guess for a single quantity based on its name
   * @param sim The simulator for which to set an inital gues
   * @param name The name of the quantity
   * @param values A C-array containing sufficient (length = number of cells) values
   */
  template < typename SimulatorT >
  void set_initial_guess(SimulatorT & sim, std::string name, double * values)
  {
    libviennashe::array_to_accessor tacc(values);
    sim.set_initial_guess(name, tacc);
  }

} // namespace libviennashe



#ifdef	__cplusplus
extern "C"
{
#endif


viennasheErrorCode viennashe_create_simulator (viennashe_simulator * sim, viennashe_device dev, viennashe_config conf)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,2,"dev");
    CHECK_ARGUMENT_FOR_NULL(conf,3,"conf");

    // Get internal configuration and device
    viennashe::config     * int_conf = reinterpret_cast<viennashe::config *>(conf);
    viennashe_device_impl * int_dev  = dev;

    // Create the internal simulator object and init with grid type and config
    viennashe_simulator_impl * int_sim  = new viennashe_simulator_impl(int_dev->device_, *int_conf);

    *sim = int_sim;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! create_simulator(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_free_simulator (viennashe_simulator_impl * sim)
{
  try
  {
    if (sim != NULL)
    {
      // The internal configuration is not destroyed!
      delete sim;
    }
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! free_simulator(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_initial_guess_from_other_sim(viennashe_simulator_impl * sim, viennashe_simulator_impl * other_sim)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(sim,1,"sim");
    CHECK_ARGUMENT_FOR_NULL(other_sim,2,"other_sim");

    //
    // Get internal structures
    viennashe_simulator_impl * int_sim       = sim;
    viennashe_simulator_impl * int_other_sim = other_sim;

    //
    // Transfer the initial guess
    libviennashe::set_initial_guess(int_sim->sim_, int_other_sim->sim_);
  }
  catch(viennashe::quantity_not_found_exception const & ex)
  {
    viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): Quantity not found exception! What? '" << ex.what() << "'" << std::endl;
    return 2;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_initial_guess(viennashe_simulator_impl * sim, const char * name, double * values)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(sim,1,"sim");

    //
    // Get internal structures
    viennashe_simulator_impl * int_sim = sim;

    //
    // Transfer the initial guess

    libviennashe::set_initial_guess(int_sim->sim_, name, values);
  }
  catch(viennashe::quantity_not_found_exception const & ex)
  {
    viennashe::log::error() << "ERROR! set_initial_guess(): Quantity not found exception! What? '" << ex.what() << "'" << std::endl;
    return 2;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! set_initial_guess(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_run(viennashe_simulator_impl * sim)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(sim,1,"sim");

    viennashe_simulator_impl * int_sim = sim;

    int_sim->sim_.run();
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! run(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


#ifdef	__cplusplus
}
#endif
