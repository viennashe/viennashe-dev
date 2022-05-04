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
/* ##################################### */
/*   This is the main header file
     for the library users               */
/* ##################################### */

#include "libviennashe/include/sys.h"
#include "libviennashe/include/error.h"

#ifndef LIBVIENNASHE_LIBVIENNASHE_H
#define	LIBVIENNASHE_LIBVIENNASHE_H
/* ##################################### */

#ifdef	__cplusplus
extern "C" {
#endif

/** @brief Initializes ViennaSHE. To be called before ViennaSHE is used */
VIENNASHE_EXPORT viennasheErrorCode viennashe_initalize(void);

/** @brief Finalizes ViennaSHE */
VIENNASHE_EXPORT viennasheErrorCode viennashe_finalize(void);

#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_LIBVIENNASHE_H */


/* ##################################### */
/*  Include all other header files here  */
/* ##################################### */

#include "libviennashe/include/material.h"
#include "libviennashe/include/device.h"
#include "libviennashe/include/config.h"
#include "libviennashe/include/simulator.h"
#include "libviennashe/include/output.h"




