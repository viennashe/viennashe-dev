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

/* C includes */
#include "libviennashe/include/sys.h"
#include "libviennashe/include/error.h"
#include "libviennashe/include/device.h"
#include "libviennashe/include/config.h"
#include "libviennashe/include/quantity.h"

#ifndef LIBVIENNASHE_SIMULATOR_H
#define	LIBVIENNASHE_SIMULATOR_H

#ifdef	__cplusplus
extern "C"
{
#endif

/*  Types  */


typedef viennashe_simulator_impl* viennashe_simulator;

/*  Functions  */

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_simulator(viennashe_simulator * sim, viennashe_device dev, viennashe_config conf);
VIENNASHE_EXPORT viennasheErrorCode viennashe_free_simulator(viennashe_simulator sim);


VIENNASHE_EXPORT viennasheErrorCode viennashe_set_initial_guess_from_other_sim(viennashe_simulator sim, viennashe_simulator other_sim);

VIENNASHE_EXPORT viennasheErrorCode viennashe_set_initial_guess(viennashe_simulator sim, const char * name, double * values);

VIENNASHE_EXPORT viennasheErrorCode viennashe_run(viennashe_simulator sim);

/*
// TODO:
//
// *) Time dependence
//
// Maybe: *) Status messages? Callbacks ??
*/


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_SIMULATOR_H */

