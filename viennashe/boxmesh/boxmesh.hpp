#ifndef VIENNASHE_BOXMESH_BOXMESH_HPP
#define VIENNASHE_BOXMESH_BOXMESH_HPP

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


#include "viennashe/forwards.h"

/** @file  viennashe/boxmesh/boxmesh.hpp
  *
  * @brief Defines a light-weight generic finite-volume datastructure.
  *
  * Decouples underlying grid from the discretization routines, so that cell-centered and vertex-centered assembly can be implemented using the same code.
  *
  */

//
// Error codes
//
#define VIENNASHE_BOXMESH_SUCCESS      0

#define VIENNASHE_BOXMESH_ERROR_GENERIC                       1
#define VIENNASHE_BOXMESH_ERROR_INVALID_BOX_OR_INTERFACE      2
#define VIENNASHE_BOXMESH_ERROR_INVALID_ELEMENT_TYPE          3
#define VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DEFINED      4
#define VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DESTROYED    5
#define VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED     6






//
// Box mesh datastructure
//

// types
typedef struct viennashe_boxmesh_impl * viennashe_boxmesh;
typedef int                             viennashe_boxmesh_element_id;

typedef enum
{
  VIENNASHE_BOXMESH_ELEMENT_INVALID = 0,
  VIENNASHE_BOXMESH_BOX,
  VIENNASHE_BOXMESH_INTERFACE
} viennashe_boxmesh_element_type;


// functions
/** @brief Creates a boxmesh and initializes its reference counter to 1.
  *
  * Note: viennashe_boxmesh_release() is called when boxmesh is destroyed.
  */
viennashe_error viennashe_boxmesh_create(viennashe_boxmesh *boxmesh);

/** @brief Manually increases the reference counter of the boxmesh by 1. */
viennashe_error viennashe_boxmesh_retain(viennashe_boxmesh boxmesh);
/** @brief Manually decreases the reference counter of the boxmesh by 1. If the reference counter reaches a value of zero, the boxmesh (including all the quantities) is destroyed */
viennashe_error viennashe_boxmesh_release(viennashe_boxmesh boxmesh);


viennashe_error viennashe_boxmesh_box_add(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id *box_id);
viennashe_error viennashe_boxmesh_interface_add(viennashe_boxmesh boxmesh,
                                                viennashe_boxmesh_element_id first_box,
                                                viennashe_boxmesh_element_id second_box,
                                                viennashe_boxmesh_element_id *interface_id);

/** @brief Returns iterators for traversing all boxes */
viennashe_error viennashe_boxmesh_boxes_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id **boxes_begin, viennashe_boxmesh_element_id **boxes_end);

/** @brief Returns iterators for traversing all interfaces in the boxmesh */
viennashe_error viennashe_boxmesh_interfaces_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id **interfaces_begin, viennashe_boxmesh_element_id **interfaces_end);

/** @brief Returns iterators for traversing all interfaces of a box */
viennashe_error viennashe_boxmesh_box_interfaces_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id box, viennashe_boxmesh_element_id **interfaces_begin, viennashe_boxmesh_element_id **interfaces_end);

/** @brief Returns iterators for traversing all boxes (one or two) of an interface */
viennashe_error viennashe_boxmesh_interface_boxes_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id interface, viennashe_boxmesh_element_id *first_box, viennashe_boxmesh_element_id *second_box);





//
// Quantities
//

// types
typedef struct viennashe_box_quantity_impl * viennashe_boxmesh_quantity;

typedef enum
{
  VIENNASHE_BOXMESH_QUANTITY_TYPE_INVALID = 0,
  VIENNASHE_BOXMESH_QUANTITY_TYPE_DOUBLE,
  VIENNASHE_BOXMESH_QUANTITY_TYPE_INT,
  VIENNASHE_BOXMESH_QUANTITY_TYPE_USER
} viennashe_boxmesh_quantity_type;

// functions

/** @brief Creates and returns a quantity on the provided boxmesh */
viennashe_error viennashe_boxmesh_quantity_create(viennashe_boxmesh boxmesh, viennashe_boxmesh_quantity *quantity, viennashe_boxmesh_element_type where_defined, const char *name, size_t bytes_per_datum);


/** @brief Returns the quantity name */
viennashe_error viennashe_boxmesh_quantity_name_get(viennashe_boxmesh_quantity quantity, const char **name);
/** @brief Returns the size of one datum */
viennashe_error viennashe_boxmesh_quantity_datum_size_get(viennashe_boxmesh_quantity quantity, size_t *bytes_per_datum);
/** @brief Returns the datum associated with the provided box or interface. Actual data type is known by the user when constructing the quantity. */
viennashe_error viennashe_boxmesh_quantity_datum_get(viennashe_boxmesh_quantity quantity, viennashe_boxmesh_element_id box_or_interface, void **datum);

/** @brief Helper routine for obtaining the quantity of a given name associated with a box mesh */
viennashe_error viennashe_boxmesh_quantity_get_by_name(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_type assoc, const char *name, viennashe_boxmesh_quantity *quantity);



#endif
