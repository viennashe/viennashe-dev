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

#include "libviennashe/include/sys.h"
#include "libviennashe/include/error.h"

#ifndef LIBVIENNASHE_CONFIG_H
#define	LIBVIENNASHE_CONFIG_H

#ifdef	__cplusplus
extern "C"
{
#endif

/*
// Types
*/

typedef struct viennashe_config_impl viennashe_config_impl; /* maps to viennashe::config */
typedef viennashe_config_impl* viennashe_config; /*! The simulator configuration */

/** @brief Enum of available linear solvers */
typedef enum { viennashe_linear_solver_dense, viennashe_linear_solver_serial,
               viennashe_linear_solver_parallel, viennashe_linear_solver_gpu_parallel } viennashe_linear_solver_id;

/** @brief Enum of available non-linear solvers*/
typedef enum { viennashe_nonlinear_solver_gummel, viennashe_nonlinear_solver_newton } viennashe_nonlinear_solver_id;

/*  Functions  */

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_config(viennashe_config * conf);
VIENNASHE_EXPORT viennasheErrorCode viennashe_free_config  (viennashe_config conf);

VIENNASHE_EXPORT viennasheErrorCode viennashe_config_standard_dd(viennashe_config conf);

VIENNASHE_EXPORT viennasheErrorCode viennashe_config_she_bipolar(viennashe_config conf, libviennashe_bool with_traps);
VIENNASHE_EXPORT viennasheErrorCode viennashe_config_she_unipolar_n(viennashe_config conf);
VIENNASHE_EXPORT viennasheErrorCode viennashe_config_she_unipolar_p(viennashe_config conf);

VIENNASHE_EXPORT viennasheErrorCode viennashe_config_enable_density_gradient(viennashe_config conf);
VIENNASHE_EXPORT viennasheErrorCode viennashe_config_disable_density_gradient(viennashe_config conf);

VIENNASHE_EXPORT viennasheErrorCode viennashe_set_linear_solver_config(viennashe_config conf, viennashe_linear_solver_id    sol_id, long max_iters);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_nonlinear_solver_config(viennashe_config conf, viennashe_nonlinear_solver_id sol_id, long max_iters, double damping);

VIENNASHE_EXPORT viennasheErrorCode viennashe_config_with_traps(viennashe_config conf, libviennashe_bool enabled);

VIENNASHE_EXPORT viennasheErrorCode viennashe_set_optical_phonon_scattering   (viennashe_config conf, libviennashe_bool enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_acoustic_phonon_scattering  (viennashe_config conf, libviennashe_bool enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_ionized_impurity_scattering (viennashe_config conf, libviennashe_bool enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_impact_ionization_scattering(viennashe_config conf, libviennashe_bool enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_electron_electron_scattering(viennashe_config conf, libviennashe_bool enabled);

/*  Getter  */

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_optical_phonon_scattering   (viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_acoustic_phonon_scattering  (viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_ionized_impurity_scattering (viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_impact_ionization_scattering(viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_electron_electron_scattering(viennashe_config conf, libviennashe_bool * enabled);

VIENNASHE_EXPORT viennasheErrorCode viennashe_config_is_with_traps(viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_config_has_density_gradient(viennashe_config conf, libviennashe_bool * enabled);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_linear_solver_config(viennashe_config conf, viennashe_linear_solver_id       * sol_id, long * max_iters);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_nonlinear_solver_config(viennashe_config conf, viennashe_nonlinear_solver_id * sol_id, long * max_iters, double * damping);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_CONFIG_H */

