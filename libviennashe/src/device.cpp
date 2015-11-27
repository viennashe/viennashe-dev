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

// C++ includes
#include "libviennashe/src/viennashe_all.hpp"
#include "viennashe/util/dump_device_mesh.hpp"
#include "viennashe/util/generate_device.hpp"

// C includes
#include "libviennashe/include/device.h"
#include "libviennashe/include/material.h"


namespace libviennashe
{

    template<typename DeviceT>
    void get_vertices_on_cell(DeviceT const & device,
                              viennashe_index_type cell_id,
                              viennashe_index_type * vertex_ids
                              )
    {
      typedef typename DeviceT::mesh_type MeshType;

      if (!vertex_ids) throw std::invalid_argument("vertices = NULL !");

      viennagrid_dimension cell_dim;
      viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

      viennagrid_element_id *cells_begin, *cells_end;
      viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

      if (/*cell_id < 0 ||*/ cell_id >= viennashe_index_type(cells_end - cells_begin)) throw std::invalid_argument("cell_id invalid !");

      viennagrid_element_id *vertices_on_cell_begin, *vertices_on_cell_end;
      viennagrid_element_boundary_elements(device.mesh(), cells_begin[cell_id], 0, &vertices_on_cell_begin, &vertices_on_cell_end);

      std::size_t j = 0;
      for (viennagrid_element_id *vit = vertices_on_cell_begin; vit != vertices_on_cell_end; ++vit, ++j)
      {
        vertex_ids[j] = viennashe_index_type(viennagrid_index_from_element_id(*vit));
      }

    } // get_vertices_on_cell


  /**
   * @brief C++ Implementation (template!) for material and doping initialization
   * @param device The ViennaSHE device
   * @param material_ids A C-array of valid ViennaSHE material ids. Has to have an entry for every cell!
   * @param doping_n A C-array of donor dopings for all cells. Cells which do not feature a doping are to have zero doping
   * @param doping_p A C-array of acceptor dopings for all cells. Cells which do not feature a doping are to have zero doping
   */
  template < typename DeviceT >
  void initalize_device(DeviceT & device, long * material_ids, double * doping_n, double * doping_p)
  {
    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

    std::size_t cell_num = cells_end - cells_begin;

    for (std::size_t i = 0; i < cell_num; ++i)
    {
      device.set_material(material_ids[i], cells_begin[i]);

      if (doping_n[i] < 0)
      {
        viennashe::log::warn() << "initalize_device(): Warning! The donor doping must be greater zero. doping_n[" << i << "] = " << doping_n[i] << std::endl;
      }
      else if (doping_p[i] < 0)
      {
        viennashe::log::warn() << "initalize_device(): Warning! The acceptor doping must be greater zero. doping_p[" << i << "] = " << doping_p[i] << std::endl;
      }
      else
      {
        // Zero doping means do not set !
        if (doping_n[i]) device.set_doping_n(doping_n[i], cells_begin[i]);
        if (doping_p[i]) device.set_doping_p(doping_p[i], cells_begin[i]);
      }
    }

  } // initalize_device

  /**
   * @brief C++ Implementation (template!) to set the contact potentials per cell
   * @param device The ViennaSHE device!
   * @param cell_ids A C-array of cell ids
   * @param values An array of contact potentials. Has to have the same length as cell_ids
   * @param len The number of cells (length of cell_ids) for which to set the contact potential
   */
  template < typename DeviceT >
  void set_contact_potential(DeviceT & device, viennashe_index_type * cell_ids, double * values, viennashe_index_type len)
  {
    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

    std::size_t cell_num = cells_end - cells_begin;

    for (std::size_t i = 0; i < len; ++i)
    {
      const viennashe_index_type id = cell_ids[i];
      if (id < cell_num)
        device.set_contact_potential(values[i], cells_begin[id] );
      else
        viennashe::log::warn() << "WARNING! set_contact_potential(): Invalid cell id '" << id << "' ... skipping!" << std::endl;
    }

  } // initalize_device

  /**
   * @brief C++ Implementation (template!) to set the contact potentials per segment
   * @param device The ViennaSHE device!
   * @param segment_id An ID to a valid ViennaGrid segment on the current mesh
   * @param value The contact potential (Volt!) to be set for all cells in the segment
   */
  template < typename DeviceT >
  void set_contact_potential(DeviceT & device, viennashe_index_type segment_id,  double value)
  {
    device.set_contact_potential(value, device.segment(static_cast<viennagrid_region_id>(segment_id)));
  } // initalize_device

} // namespace libviennashe

#ifdef	__cplusplus
extern "C" {
#endif


viennasheErrorCode viennashe_create_1d_device(viennashe_device * dev, double len_x, size_t points_x)
{
  try
  {
    viennashe_device_impl * int_dev = new viennashe_device_impl();

    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev, 1, "dev");
    if (len_x <= 0.0)
    {
      viennashe::log::error() << "ERROR! viennashe_create_1d_device(): len_x must be greater zero. len_x = " << len_x << std::endl;
      return 2;
    }
    if (points_x == 0)
    {
      viennashe::log::error() << "ERROR! viennashe_create_1d_device(): points_x must be greater zero. points_x = " << points_x << std::endl;
      return 3;
    }


    viennashe::util::device_generation_config generator_params;
    generator_params.add_segment(0.0, len_x, static_cast<unsigned long>(points_x));
    int_dev->device_.generate_mesh(generator_params);

    *dev = int_dev;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_1d_device(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


viennasheErrorCode viennashe_free_device(viennashe_device dev)
{
  try
  {
    if (dev)
      delete dev;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! free_device(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_initalize_device(viennashe_device dev, viennashe_material_id * material_ids, double * doping_n, double * doping_p)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    if (doping_n == 0 || doping_p == 0 || material_ids == 0)
    {
      viennashe::log::error() << "ERROR! initalize_device(): The arrays material_ids, doping_n and doping_p must be arrays " << std::endl;
      return 2;
    }

    viennashe_device_impl * int_dev = dev;

    libviennashe::initalize_device(int_dev->device_, material_ids, doping_n, doping_p);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! initalize_device(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_material_on_segment(viennashe_device dev, viennashe_material_id material_id, viennashe_index_type segment_id)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    if (material_id <  0)
    {
      viennashe::log::error() << "ERROR! viennashe_set_material_on_segment(): The material_id needs to be greater zero. material_id = " << material_id  << std::endl;
      return 2;
    }

    viennashe_device_impl * int_dev = dev;

    viennagrid_region region;
    viennagrid_mesh_region_get(int_dev->device_.mesh(), segment_id, &region);
    int_dev->device_.set_material(material_id, region);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_material_on_segment(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_doping_n_on_segment(viennashe_device dev, double doping_n, viennashe_index_type segment_id )
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    if (doping_n <  0)
    {
      viennashe::log::error() << "ERROR! viennashe_set_doping_n_on_segment(): The doping_n needs to be greater zero. doping_n = " << doping_n  << std::endl;
      return 2;
    }

    viennashe_device_impl * int_dev = dev;

    viennagrid_region region;
    viennagrid_mesh_region_get(int_dev->device_.mesh(), segment_id, &region);
    int_dev->device_.set_doping_n(doping_n, region);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_doping_n_on_segment(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_doping_p_on_segment(viennashe_device dev, double doping_p, viennashe_index_type segment_id )
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    if (doping_p <  0)
    {
      viennashe::log::error() << "ERROR! viennashe_set_doping_p_on_segment(): The doping_p needs to be greater zero. doping_p = " << doping_p  << std::endl;
      return 2;
    }

    viennashe_device_impl * int_dev = dev;

    viennagrid_region region;
    viennagrid_mesh_region_get(int_dev->device_.mesh(), segment_id, &region);
    int_dev->device_.set_doping_p(doping_p, region);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_doping_p_on_segment(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_contact_potential_cells(viennashe_device dev, viennashe_index_type * cell_ids, double * values, viennashe_index_type len)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(cell_ids,2,"cell_ids");
    CHECK_ARGUMENT_FOR_NULL(values,3,"values");

    viennashe_device_impl * int_dev = dev;

    libviennashe::set_contact_potential(int_dev->device_, cell_ids, values, len);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_set_contact_potential_cells(): UNKOWN ERROR!" << std::endl;
    return 1;
  }
  return 0;
}

viennasheErrorCode viennashe_set_contact_potential_segment(viennashe_device dev, double value, viennashe_index_type segment_id)
{
  try
  {
    //
    // Checks
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");

    viennashe_device_impl * int_dev = dev;

    libviennashe::set_contact_potential(int_dev->device_, segment_id, value);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! set_contact_potential(): UNKOWN ERROR!" << std::endl;
    return 1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_vertices(viennashe_device dev, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,2,"num");

    viennashe_device_impl * int_dev = dev;

    viennagrid_int vertex_count;
    viennagrid_mesh_element_count(int_dev->device_.mesh(), 0, &vertex_count);
    *num = vertex_count;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! get_num_vertices(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_cells(viennashe_device dev, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,2,"num");

    viennashe_device_impl * int_dev = dev;

    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(int_dev->device_.mesh(), &cell_dim);

    viennagrid_int cell_count;
    viennagrid_mesh_element_count(int_dev->device_.mesh(), cell_dim, &cell_count);
    *num = cell_count;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! get_num_cells(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_segments(viennashe_device dev, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,2,"num");

    viennashe_device_impl * int_dev = dev;

    viennagrid_region_id *regions_begin, *regions_end;
    viennagrid_mesh_regions_get(int_dev->device_.mesh(), &regions_begin, &regions_end);
    *num = regions_end - regions_begin;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! get_num_segments(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_vertices_on_segment (viennashe_device dev, viennashe_index_type segment_id, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,3,"num");

    //
    viennashe_device_impl * int_dev = dev;

    viennagrid_element_id *vertices_begin, *vertices_end;
    viennagrid_mesh_elements_get(int_dev->device_.mesh(), 0, &vertices_begin, &vertices_end);

    viennagrid_region region = int_dev->device_.segment(segment_id);

    viennashe_index_type vertex_num = 0;

    for (viennagrid_element_id *vit = vertices_begin; vit != vertices_end; ++vit)
    {
      viennagrid_bool vertex_in_region;
      viennagrid_region_contains_element(region, *vit, &vertex_in_region);
      if (vertex_in_region)
        ++vertex_num;
    }

    *num = vertex_num;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! get_num_vertices_on_segment(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_cells_on_segment(viennashe_device dev, viennashe_index_type segment_id, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,3,"num");

    viennashe_device_impl * int_dev = dev;

    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(int_dev->device_.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(int_dev->device_.mesh(), cell_dim, &cells_begin, &cells_end);

    viennagrid_region region = int_dev->device_.segment(segment_id);

    viennashe_index_type cell_num = 0;

    for (viennagrid_element_id *cit = cells_begin; cit != cells_end; ++cit)
    {
      viennagrid_bool cell_in_region;
      viennagrid_region_contains_element(region, *cit, &cell_in_region);
      if (cell_in_region)
        ++cell_num;
    }

    *num = cell_num;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! get_num_cells_on_segment(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_dimension (viennashe_device dev, viennashe_index_type * dim)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(dim,2,"dim");

    viennashe_device_impl * int_dev = dev;

    viennagrid_dimension geo_dim;
    viennagrid_mesh_geometric_dimension_get(int_dev->device_.mesh(), &geo_dim);
    *dim = geo_dim;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_dimension(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_num_vertices_per_cell (viennashe_device dev, viennashe_index_type * num)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(num,2,"num");

    viennashe_device_impl * int_dev = dev;

    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(int_dev->device_.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(int_dev->device_.mesh(), cell_dim, &cells_begin, &cells_end);

    viennagrid_element_id *vertices_on_cell_begin, *vertices_on_cell_end;
    viennagrid_element_boundary_elements(int_dev->device_.mesh(), cells_begin[0], 0, &vertices_on_cell_begin, &vertices_on_cell_end);
    *num = vertices_on_cell_end - vertices_on_cell_begin;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_num_vertices_per_cell(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device(viennashe_device * dev, viennashe_index_type geo_dim,
                                                            double ** vertices, viennashe_index_type num_vertices,
                                                            viennashe_index_type ** cells, viennashe_index_type num_cells,
                                                            viennashe_index_type * segmentation)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(vertices,3,"vertices");
    CHECK_ARGUMENT_FOR_NULL(cells,5,"cells");

    viennashe_device_impl * int_dev = new viennashe_device_impl();

    viennashe::util::device_from_array_generator<viennashe_index_type> gen(geo_dim, vertices, cells, segmentation, num_vertices, num_cells);
    int_dev->device_.generate_mesh(gen);

    *dev = int_dev;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}



VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_flat(viennashe_device * dev, viennashe_index_type geo_dim,
                                                                 double * vertices, viennashe_index_type num_vertices,
                                                                 viennashe_index_type * cells, viennashe_index_type num_cells,
                                                                 viennashe_index_type * segmentation)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(vertices,3,"vertices");
    CHECK_ARGUMENT_FOR_NULL(cells,5,"cells");

    viennashe_device_impl * int_dev = new viennashe_device_impl();

    viennashe::util::device_from_flat_array_generator<viennashe_index_type> gen(geo_dim, vertices, cells, segmentation, num_vertices, num_cells);
    int_dev->device_.generate_mesh(gen);

    *dev = int_dev;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device_flat(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_from_file(viennashe_device * dev, char const * filename)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(filename,2,"filename");

    viennashe_device_impl * int_dev = new viennashe_device_impl();
    int_dev->device_.load_mesh(filename);

    *dev = int_dev;
  }
  catch(std::exception const & ex)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device_from_file(): " << ex.what() << std::endl;
    return -1;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device_from_file(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_grid (viennashe_device  dev, double ** vertices, viennashe_index_type * num_vertices,
                                       viennashe_index_type ** cells, viennashe_index_type * num_cells)
{
  try
  {
    //
    // CHECKS
    if (dev == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_get_grid(): Mesh must exist (dev != NULL) "<< std::endl;
      return 1;
    }
    if (num_vertices == NULL || vertices == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_get_grid(): num_vertices must be valid. num_vertices = "
                              << num_vertices << " and vertices = " << vertices << std::endl;
      return 2;
    }
    if (num_cells == NULL || cells == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_get_grid(): num_cells must be valid. num_cells = "
                              << num_cells << " and cells = " << cells << std::endl;
      return 4;
    }

    viennashe::util::dump_mesh(dev->device_, vertices, *num_vertices, cells, *num_cells);
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_grid(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
} // viennashe_get_grid



viennasheErrorCode viennashe_get_nth_vertex(viennashe_device dev, viennashe_index_type vid, double * x, double * y, double * z)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(x,3,"x");
    CHECK_ARGUMENT_FOR_NULL(y,4,"y");
    CHECK_ARGUMENT_FOR_NULL(z,5,"z");

    viennashe_index_type num_vertices = 0;
    viennashe_get_num_vertices(dev, &num_vertices);

    if ( static_cast<long>(vid) < 0 || vid >= num_vertices)  // static_cast to silence tautology warnings (vid < 0 always fulfilled for unsigned integers)
    {
      viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): The given index id (vid) must not be smaller 0 or greater " << num_vertices << std::endl;
      return 2;
    }

    *x = 0; *y = 0; *z = 0;

    viennagrid_element_id *vertices_begin, *vertices_end;
    viennagrid_mesh_elements_get(dev->device_.mesh(), 0, &vertices_begin, &vertices_end);

    viennagrid_numeric *coords;
    viennagrid_mesh_vertex_coords_get(dev->device_.mesh(), vertices_begin[vid], &coords);

    viennagrid_dimension geo_dim;
    viennagrid_mesh_geometric_dimension_get(dev->device_.mesh(), &geo_dim);
    *x = coords[0];
    if (geo_dim > 1)
      *y = coords[1];
    if (geo_dim > 2)
      *z = coords[2];

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_nth_cell  (viennashe_device dev, viennashe_index_type cid, viennashe_index_type * vertex_id_list)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(vertex_id_list,3,"vertex_id_list");

    viennashe_index_type num_cells = 0;
    viennashe_get_num_cells(dev, &num_cells);

    if (static_cast<long>(cid) < 0 || cid >= num_cells)  // static_cast to silence tautology warnings (vid < 0 always fulfilled for unsigned integers)
    {
      viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): The given index id (cid) must not be smaller 0 or greater " << num_cells << std::endl;
      return 2;
    }

    libviennashe::get_vertices_on_cell(dev->device_, cid, vertex_id_list);

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


#ifdef	__cplusplus
}
#endif
