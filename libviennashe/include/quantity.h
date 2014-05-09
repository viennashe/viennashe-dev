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

#ifndef LIBVIENNASHE_QUANTITY_H
#define	LIBVIENNASHE_QUANTITY_H

/* C includes */
#include "libviennashe/include/sys.h"
#include "libviennashe/include/error.h"
#include "libviennashe/include/simulator.h"

#ifdef	__cplusplus
extern "C"
{
#endif

/* ******************************************* */
/*           Type definitions                  */

typedef viennashe_quan_register_impl* viennashe_quan_register; /*! The quantity registry */

/** @brief Enum of available charge carrier types */
typedef enum { viennashe_electron_id, viennashe_hole_id } viennashe_carrier_ids;

/* ******************************************* */


VIENNASHE_EXPORT viennasheErrorCode viennashe_create_quantity_register(viennashe_quan_register * reg, viennashe_simulator sim);
VIENNASHE_EXPORT viennasheErrorCode viennashe_free_quantity_register(viennashe_quan_register reg);

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_cell_based(viennashe_quan_register reg, viennashe_index_type * num);

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_cell_based_quantity_list(viennashe_quan_register reg, char ** names);

VIENNASHE_EXPORT viennasheErrorCode viennashe_has_cell_based_quantity(viennashe_quan_register reg, char const * name, libviennashe_bool * exists);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_cell_based_quantity(viennashe_quan_register reg, char const * name, double ** values, viennashe_index_type * len);


/* ************** */
/* SHE quantities */
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_edf            (viennashe_quan_register reg, viennashe_carrier_ids ctype, double ** energies, double ** values, viennashe_index_type * len);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_dos            (viennashe_quan_register reg, viennashe_carrier_ids ctype, double ** energies, double ** values, viennashe_index_type * len);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_group_velocity (viennashe_quan_register reg, viennashe_carrier_ids ctype, double ** energies, double ** values, viennashe_index_type * len);
/* ************** */


/* ******************************************* */
/* Convenience functions for memory management */

VIENNASHE_EXPORT viennasheErrorCode viennashe_prealloc_cell_based_quantity(viennashe_device dev, double *** uarray, viennashe_index_type ** len);

VIENNASHE_EXPORT viennasheErrorCode viennashe_free_cell_based_quantity(viennashe_device dev, double *** uarray, viennashe_index_type ** len);

/* ******************************************* */


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_QUANTITY_H */

