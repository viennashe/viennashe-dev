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

/* Internally used for correct export statements */
#ifdef WIN32 /* We need the dllexport statement for Windows only */
#define VIENNASHE_EXPORT __declspec(dllexport)
#else /* No export statement for Unix-based systems */
#define VIENNASHE_EXPORT
#endif

#ifndef LIBVIENNASHE_SYS_H
#define	LIBVIENNASHE_SYS_H

#ifdef	__cplusplus
extern "C"
{
#endif

/* @brief Implementation of the C++ bool variable for C98 and C90 */
typedef enum { libviennashe_false=0, libviennashe_true } libviennashe_bool;

/** @brief Device implementation type */
typedef struct viennashe_device_impl     viennashe_device_impl;

/** @brief Simulator implementation type */
typedef struct viennashe_simulator_impl  viennashe_simulator_impl;

/** @brief Quantity register implementation type */
typedef struct viennashe_quan_register_impl viennashe_quan_register_impl;

typedef unsigned long  viennashe_index_type;

#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_SYS_H */

