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
    viennashe::config                * int_conf = reinterpret_cast<viennashe::config *>(conf);
    viennashe_device_impl * int_dev  = dev;

    // Check the device
    if (!int_dev->is_valid())
    {
      viennashe::log::error() << "ERROR! create_simulator(): The device (dev) must be valid!" << std::endl;
      return 2;
    }
    // Create the internal simulator object and init with grid type and config
    viennashe_simulator_impl * int_sim  = new viennashe_simulator_impl(int_dev->stype, int_conf);

    //
    // Create viennashe::simulator per grid type
    if(int_dev->stype == libviennashe::meshtype::line_1d)
    {
      int_sim->sim1d  = new viennashe_simulator_impl::sim1d_type(*(int_dev->device_1d), *int_conf);
    }
    else if(int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      int_sim->simq2d  = new viennashe_simulator_impl::simq2d_type(*(int_dev->device_quad_2d), *int_conf);
    }
    else if(int_dev->stype == libviennashe::meshtype::triangular_2d)
    {
      int_sim->simt2d  = new viennashe_simulator_impl::simt2d_type(*(int_dev->device_tri_2d), *int_conf);
    }
    else if(int_dev->stype == libviennashe::meshtype::hexahedral_3d)
    {
      int_sim->simh3d  = new viennashe_simulator_impl::simh3d_type(*(int_dev->device_hex_3d), *int_conf);
    }
    else if(int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      int_sim->simt3d  = new viennashe_simulator_impl::simt3d_type(*(int_dev->device_tet_3d), *int_conf);
    }
    else
    {
      viennashe::log::error() << "ERROR! create_device(): The given mesh is malconfigured!" << std::endl;
      return 2;
    }
    // RETURN
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
    viennashe_simulator_impl * int_sim       = (sim);
    viennashe_simulator_impl * int_other_sim = (other_sim);
    //
    // Checks
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }
    if (!int_other_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): The simulator (other_sim) must be valid!" << std::endl;
      return 1;
    }
    if (int_sim->stype != int_other_sim->stype)
    {
      viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): The simulators (sim and other_sim)"
                              << " must be operating on the same grid!" << std::endl;
      return 2;
    }

    //
    // Transfer the initial guess

    if(int_sim->stype == libviennashe::meshtype::line_1d)
    {
      libviennashe::set_initial_guess(*(int_sim->sim1d), *(int_other_sim->sim1d));
    }
    else if(int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      libviennashe::set_initial_guess(*(int_sim->simq2d), *(int_other_sim->simq2d));
    }
    else if(int_sim->stype == libviennashe::meshtype::triangular_2d)
    {
      libviennashe::set_initial_guess(*(int_sim->simt2d), *(int_other_sim->simt2d));
    }
    else if(int_sim->stype == libviennashe::meshtype::hexahedral_3d)
    {
      libviennashe::set_initial_guess(*(int_sim->simh3d), *(int_other_sim->simh3d));
    }
    else if(int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      libviennashe::set_initial_guess(*(int_sim->simt3d), *(int_other_sim->simt3d));
    }
    else
    {
      viennashe::log::error() << "ERROR! set_initial_guess_from_other_sim(): Unkown grid type!" << std::endl;
      return 1;
    }
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
    // Checks
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! set_initial_guess(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    //
    // Transfer the initial guess

    if(int_sim->stype == libviennashe::meshtype::line_1d)
    {
      libviennashe::set_initial_guess(*(int_sim->sim1d), name, values);
    }
    else if(int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      libviennashe::set_initial_guess(*(int_sim->simq2d), name, values);
    }
    else if(int_sim->stype == libviennashe::meshtype::triangular_2d)
    {
      libviennashe::set_initial_guess(*(int_sim->simt2d), name, values);
    }
    else if(int_sim->stype == libviennashe::meshtype::hexahedral_3d)
    {
      libviennashe::set_initial_guess(*(int_sim->simh3d), name, values);
    }
    else if(int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      libviennashe::set_initial_guess(*(int_sim->simt3d), name, values);
    }
    else
    {
      viennashe::log::error() << "ERROR! set_initial_guess(): Unkown grid type!" << std::endl;
      return -2;
    }
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

    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! run(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    if(int_sim->stype == libviennashe::meshtype::line_1d)
    {
      int_sim->sim1d->run();
    }
    else if(int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      int_sim->simq2d->run();
    }
    else if(int_sim->stype == libviennashe::meshtype::triangular_2d)
    {
      int_sim->simt2d->run();
    }
    else if(int_sim->stype == libviennashe::meshtype::hexahedral_3d)
    {
      int_sim->simh3d->run();
    }
    else if(int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      int_sim->simt3d->run();
    }
    else
    {
      viennashe::log::error() << "ERROR! run(): Unkown grid type!" << std::endl;
      return -2;
    }
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
