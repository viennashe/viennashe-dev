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
#include "libviennashe/include/material.h"

#ifndef LIBVIENNASHE_DEVICE_H
#define	LIBVIENNASHE_DEVICE_H

#ifdef	__cplusplus
extern "C" {
#endif

/*
// Types
*/

/** @brief Enum of available toplogical mesh configurations */
typedef enum { viennashe_line_1d,
               viennashe_quadrilateral_2d, viennashe_triangular_2d,
               viennashe_hexahedral_3d, viennashe_tetrahedral_3d }
        viennashe_topology_type_id;


typedef viennashe_device_impl*        viennashe_device; /*! The device! */

/*
// Functions
*/

VIENNASHE_EXPORT viennasheErrorCode viennashe_free_device(viennashe_device dev);

VIENNASHE_EXPORT viennasheErrorCode viennashe_initalize_device(viennashe_device dev, viennashe_material_id * material_ids, double * doping_n, double * doping_p);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_material_on_segment(viennashe_device dev, viennashe_material_id material_id, viennashe_index_type segment_id);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_doping_n_on_segment(viennashe_device dev, double doping_n, viennashe_index_type segment_id);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_doping_p_on_segment(viennashe_device dev, double doping_p, viennashe_index_type segment_id);

VIENNASHE_EXPORT viennasheErrorCode viennashe_set_contact_potential_cells(viennashe_device dev, viennashe_index_type * cell_ids, double * values, viennashe_index_type len);
VIENNASHE_EXPORT viennasheErrorCode viennashe_set_contact_potential_segment(viennashe_device dev, double   value, viennashe_index_type   segment_id);


/* *************** */
/* Device creators */
/* *************** */

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_1d_device(viennashe_device * dev, double len_x, size_t points_x);

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device(viennashe_device * dev, viennashe_topology_type_id topology_id,
                                                            double ** vertices, viennashe_index_type num_vertices,
                                                            viennashe_index_type ** cells, viennashe_index_type num_cells,
                                                            viennashe_index_type * segmentation);

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_flat(viennashe_device * dev, viennashe_topology_type_id topology_id,
                                                                 double * vertices, viennashe_index_type num_vertices,
                                                                 viennashe_index_type * cells, viennashe_index_type num_cells,
                                                                 viennashe_index_type * segmentation);

VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_from_file(viennashe_device * dev, viennashe_topology_type_id topology_id, char const * filename);


/* *************** */
/*   Mesh getter   */
/* *************** */

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_vertices(viennashe_device dev, viennashe_index_type * num);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_cells   (viennashe_device dev, viennashe_index_type * num);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_segments(viennashe_device dev, viennashe_index_type * num);

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_vertices_on_segment(viennashe_device dev, viennashe_index_type segment_id, viennashe_index_type * num);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_cells_on_segment   (viennashe_device dev, viennashe_index_type segment_id, viennashe_index_type * num);

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_dimension            (viennashe_device dev, viennashe_index_type * dim);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_num_vertices_per_cell(viennashe_device dev, viennashe_index_type * num);

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_grid(viennashe_device dev, double ** vertices, viennashe_index_type * num_vertices,
                                                       viennashe_index_type ** cells, viennashe_index_type * num_cells);

/* ****************** */
/*   Mesh iterators   */
/* ****************** */

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_nth_vertex(viennashe_device dev, viennashe_index_type vid, double * x, double * y, double * z);
VIENNASHE_EXPORT viennasheErrorCode viennashe_get_nth_cell  (viennashe_device dev, viennashe_index_type cid, viennashe_index_type * vertex_id_list);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBVIENNASHE_DEVICE_H */

