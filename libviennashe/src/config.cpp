/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------
 * s
                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

// C++ includes
#include "viennashe_all.hpp"

// C includes
#include "libviennashe/include/config.h"


#ifndef __FUNCTION_NAME__
    #define LIBVIENNASHE_FUNCTION_NAME_MARKER
    #ifdef WIN32   // WINDOWS
        #define __FUNCTION_NAME__   __FUNCTION__
    #else          // UNIX
        #define __FUNCTION_NAME__   __func__
    #endif
#endif

/* DEFINE THE CONFIG GETTER - Undefined at the end of this file */
#define LIBVIENNASHE_CONFIG_GETTER(conf, retval, func_name, code) try { if ((conf) == NULL) \
                                           { \
                                           viennashe::log::error() << "ERROR! " << (func_name) << "(): config must not be NULL! " << std::endl; return 1; \
                                           } \
                                           if ((retval) == NULL) \
                                           { \
                                           viennashe::log::error() << "ERROR! " << (func_name) << "(): The return pointer must not be NULL! " << std::endl; return 2; \
                                           } \
                                           viennashe::config * int_conf = reinterpret_cast<viennashe::config *>((conf)); \
                                           code; \
                                           } catch(...) {\
                                           viennashe::log::error() << "ERROR! " << (func_name) << "(): UNKOWN ERROR!" << std::endl;  return -1;\
                                           } return 0;

#ifdef	__cplusplus
extern "C"
{
#endif

/**
 *
 * @return
 */
viennasheErrorCode viennashe_create_config (viennashe_config * conf)
{
  try
  {
    // Create the config
    viennashe::config * int_conf = new viennashe::config();
    // RETURN
    *conf = reinterpret_cast<viennashe_config>(int_conf);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! create_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_free_config (viennashe_config conf)
{
  try
  {
    if (conf != NULL)
    {
      delete reinterpret_cast<viennashe::config *>(conf);
    }
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! free_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_config_standard_dd (viennashe_config conf)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);

    int_conf->with_electrons(true);
    int_conf->with_holes(true);
    int_conf->with_traps(false);
    int_conf->set_hole_equation(viennashe::EQUATION_CONTINUITY);
    int_conf->set_electron_equation(viennashe::EQUATION_CONTINUITY);
    int_conf->quantum_correction(false);
    int_conf->with_quantum_correction(false);

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_standard_dd(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @param with_traps
 * @return
 */
viennasheErrorCode viennashe_config_she_bipolar (viennashe_config conf, libviennashe_bool with_traps)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->with_electrons(true);
    int_conf->with_holes(true);

    if (with_traps) int_conf->with_traps(true);
    else            int_conf->with_traps(false);

    int_conf->set_hole_equation(viennashe::EQUATION_SHE);
    int_conf->set_electron_equation(viennashe::EQUATION_SHE);
    //int_conf->quantum_correction(false);
    //int_conf->with_quantum_correction(false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_she_bipolar(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_config_she_unipolar_n (viennashe_config conf )
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->with_electrons(true);
    int_conf->with_holes(true);

    int_conf->with_traps(false);

    int_conf->set_hole_equation(viennashe::EQUATION_CONTINUITY);
    int_conf->set_electron_equation(viennashe::EQUATION_SHE);
    //int_conf->quantum_correction(false);
    //int_conf->with_quantum_correction(false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_she_unipolar_n(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_config_she_unipolar_p (viennashe_config conf )
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->with_electrons(true);
    int_conf->with_holes(true);

    int_conf->with_traps(false);

    int_conf->set_electron_equation(viennashe::EQUATION_CONTINUITY);
    int_conf->set_hole_equation(viennashe::EQUATION_SHE);
    //int_conf->quantum_correction(false);
    //int_conf->with_quantum_correction(false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_she_unipolar_p(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_config_enable_density_gradient(viennashe_config conf)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->quantum_correction(true);
    int_conf->with_quantum_correction(true);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_enable_density_gradient(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @return
 */
viennasheErrorCode viennashe_config_disable_density_gradient(viennashe_config conf)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->quantum_correction(false);
    int_conf->with_quantum_correction(false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! config_enable_density_gradient(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @param sol_id
 * @param max_iters
 * @return
 */
viennasheErrorCode viennashe_set_linear_solver_config    (viennashe_config conf, viennashe_linear_solver_id    sol_id, long max_iters )
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    if (max_iters <= 0)
    {
      viennashe::log::error() << "ERROR! set_linear_solver_config(): max_iters must be greater zero. max_iters = " << max_iters << std::endl;
      return 3;
    }

    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    switch(sol_id)
    {
      case viennashe_linear_solver_dense:
        int_conf->linear_solver().set(viennashe::config::linear_solver_config_type::dense_linear_solver);
        break;
      case viennashe_linear_solver_serial:
        int_conf->linear_solver().set(viennashe::config::linear_solver_config_type::serial_linear_solver);
        break;
      case viennashe_linear_solver_parallel:
        int_conf->linear_solver().set(viennashe::config::linear_solver_config_type::parallel_linear_solver);
        break;
      case viennashe_linear_solver_gpu_parallel:
        int_conf->linear_solver().set(viennashe::config::linear_solver_config_type::gpu_parallel_linear_solver);
        break;
      default:
        viennashe::log::error() << "ERROR! set_linear_solver_config(): sol_id must be a valid solver id!" << std::endl;
        return 2;
    }
    int_conf->linear_solver().max_iters(static_cast<std::size_t>(max_iters));
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! set_linear_solver_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

/**
 *
 * @param conf
 * @param sol_id
 * @param max_iters
 * @param damping
 * @return
 */
viennasheErrorCode viennashe_set_nonlinear_solver_config (viennashe_config conf, viennashe_nonlinear_solver_id sol_id, long max_iters, double damping)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    if (max_iters <= 0)
    {
      viennashe::log::error() << "ERROR! set_nonlinear_solver_config(): max_iters must be greater zero. max_iters = " << max_iters << std::endl;
      return 3;
    }
    if (damping <= 0 || damping > 1)
    {
      viennashe::log::error() << "ERROR! set_nonlinear_solver_config(): damping must be in (0,1]. damping = " << damping << std::endl;
      return 4;
    }

    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    switch(sol_id)
    {
      case viennashe_nonlinear_solver_gummel:
        int_conf->nonlinear_solver().set(viennashe::config::nonlinear_solver_config_type::gummel_nonlinear_solver);
        break;
      case viennashe_nonlinear_solver_newton:
        int_conf->nonlinear_solver().set(viennashe::config::nonlinear_solver_config_type::newton_nonlinear_solver);
        break;
      default:
        viennashe::log::error() << "ERROR! set_nonlinear_solver_config(): sol_id must be a valid solver id!" << std::endl;
        return 2;
    }
    int_conf->nonlinear_solver().max_iters(static_cast<std::size_t>(max_iters));
    int_conf->nonlinear_solver().damping(damping);

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! set_nonlinear_solver_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


viennasheErrorCode viennashe_config_with_traps(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->with_traps((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_config_with_traps(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


viennasheErrorCode viennashe_set_optical_phonon_scattering(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->scattering().optical_phonon().enabled((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_optical_phonon_scattering(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_acoustic_phonon_scattering(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->scattering().acoustic_phonon().enabled((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_acoustic_phonon_scattering(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_ionized_impurity_scattering(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->scattering().ionized_impurity().enabled((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_ionized_impurity_scattering(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_impact_ionization_scattering(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->scattering().impact_ionization().enabled((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_impact_ionization_scattering(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_electron_electron_scattering(viennashe_config conf, libviennashe_bool enabled)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    viennashe::config * int_conf = reinterpret_cast<viennashe::config *>(conf);
    int_conf->scattering().electron_electron((enabled != 0) ? true : false);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_electron_electron_scattering(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}



/* +++++++++++++++ GETTER +++++++++++++++ */

viennasheErrorCode viennashe_get_optical_phonon_scattering   (viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->scattering().optical_phonon().enabled() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_get_acoustic_phonon_scattering  (viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->scattering().acoustic_phonon().enabled() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_get_ionized_impurity_scattering (viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->scattering().ionized_impurity().enabled() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_get_impact_ionization_scattering(viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->scattering().impact_ionization().enabled() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_get_electron_electron_scattering(viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->scattering().electron_electron() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_config_is_with_traps(viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = int_conf->with_traps() ? libviennashe_true : libviennashe_false)
}

viennasheErrorCode viennashe_config_has_density_gradient(viennashe_config conf, libviennashe_bool * enabled)
{
LIBVIENNASHE_CONFIG_GETTER(conf, enabled, __FUNCTION_NAME__, *enabled = ((int_conf->quantum_correction() && int_conf->with_quantum_correction()) ? libviennashe_true : libviennashe_false) )
}

viennasheErrorCode viennashe_get_linear_solver_config(viennashe_config conf, viennashe_linear_solver_id       * sol_id, long * max_iters)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    CHECK_ARGUMENT_FOR_NULL(sol_id,2,"sol_id");
    CHECK_ARGUMENT_FOR_NULL(max_iters,3,"max_iters");

    viennashe::config * int_conf = reinterpret_cast<viennashe::config *> ((conf));
    if (int_conf->linear_solver().id() == viennashe::config::linear_solver_config_type::dense_linear_solver) *sol_id = viennashe_linear_solver_dense;
    else if (int_conf->linear_solver().id() == viennashe::config::linear_solver_config_type::serial_linear_solver) *sol_id = viennashe_linear_solver_serial;
    else if (int_conf->linear_solver().id() == viennashe::config::linear_solver_config_type::parallel_linear_solver) *sol_id = viennashe_linear_solver_parallel;
    else if (int_conf->linear_solver().id() == viennashe::config::linear_solver_config_type::gpu_parallel_linear_solver) *sol_id = viennashe_linear_solver_gpu_parallel;

    *max_iters = static_cast<long>(int_conf->linear_solver().max_iters());

  }
  catch (...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_linear_solver_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_nonlinear_solver_config(viennashe_config conf, viennashe_nonlinear_solver_id * sol_id, long * max_iters, double * damping)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(conf,1,"conf");
    CHECK_ARGUMENT_FOR_NULL(sol_id,2,"sol_id");
    CHECK_ARGUMENT_FOR_NULL(max_iters,3,"max_iters");
    CHECK_ARGUMENT_FOR_NULL(damping,4,"damping");

    viennashe::config * int_conf = reinterpret_cast<viennashe::config *> ((conf));
    if (int_conf->nonlinear_solver().id() == viennashe::config::nonlinear_solver_config_type::gummel_nonlinear_solver) *sol_id = viennashe_nonlinear_solver_gummel;
    else if (int_conf->nonlinear_solver().id() == viennashe::config::nonlinear_solver_config_type::newton_nonlinear_solver) *sol_id = viennashe_nonlinear_solver_newton;

    *max_iters = static_cast<long>(int_conf->nonlinear_solver().max_iters());
    *damping = int_conf->nonlinear_solver().damping();

  }
  catch (...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_nonlinear_solver_config(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


#ifdef	__cplusplus
}
#endif

/* UNDEFINE THE CONFIG GETTER */
#undef LIBVIENNASHE_CONFIG_GETTER

#ifdef LIBVIENNASHE_FUNCTION_NAME_MARKER
#undef __FUNCTION_NAME__
#undef LIBVIENNASHE_FUNCTION_NAME_MARKER
#endif

