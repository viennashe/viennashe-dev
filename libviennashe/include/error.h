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

#include "libviennashe/include/sys.h"

#ifndef  LIBVIENNASHE_ERROR_H
#define	 LIBVIENNASHE_ERROR_H

#ifdef	__cplusplus
extern "C"
{
#endif

typedef int viennasheErrorCode; /*! The error code type */

/** @brief Error code interpretation function */
VIENNASHE_EXPORT viennasheErrorCode viennashe_error(viennasheErrorCode ecode);

#define LIBVIENNASHE_CHKERR(x)  if((x)) return viennashe_error((x));
#define LIBVIENNASHE_NO_ERROR 0

#ifdef	__cplusplus
}
#endif

#endif	/*  LIBVIENNASHE_ERROR_H */

