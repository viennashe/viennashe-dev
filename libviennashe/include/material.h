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

#ifndef LIBVIENNASHE_MATERIAL_H
#define	LIBVIENNASHE_MATERIAL_H

#ifdef	__cplusplus
extern "C"
{
#endif

typedef long viennashe_material_id; /*! The material id type */

/*   Material related stuff   */


VIENNASHE_EXPORT viennasheErrorCode viennashe_get_silicon_id(viennashe_material_id * id);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_metal_id  (viennashe_material_id * id);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_sio2_id   (viennashe_material_id * id);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_hfo2_id   (viennashe_material_id * id);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_MATERIAL_H */

