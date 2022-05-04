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

#ifndef LIBVIENNASHE_OUTPUT_H
#define	LIBVIENNASHE_OUTPUT_H

/* C includes */
#include "libviennashe/include/sys.h"
#include "libviennashe/include/device.h"
#include "libviennashe/include/config.h"
#include "libviennashe/include/simulator.h"
#include "libviennashe/include/quantity.h"

#ifdef	__cplusplus
extern "C"
{
#endif

VIENNASHE_EXPORT  viennasheErrorCode viennashe_write_to_gnuplot(viennashe_quan_register reg, char const * name, char const * filename);

VIENNASHE_EXPORT  viennasheErrorCode viennashe_write_she_results_to_gnuplot(viennashe_quan_register reg, viennashe_carrier_ids ctype, char const * filename);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_OUTPUT_H */

