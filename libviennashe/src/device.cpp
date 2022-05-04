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

// C++ includes
#include "libviennashe/src/viennashe_all.hpp"
#include "viennashe/util/dump_device_mesh.hpp"
#include "viennashe/util/generate_device.hpp"

// C includes
#include "libviennashe/include/device.h"
#include "libviennashe/include/material.h"


namespace libviennashe
{

    template < typename DeviceT >
    void get_vertices_on_cell(DeviceT const & device,
                              viennashe_index_type cell_id,
                              viennashe_index_type * vertex_ids
                              )
    {
      typedef typename DeviceT::mesh_type MeshType;
      typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;
      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::const_vertex_range<CellType>::type     VertexOnCellContainer;
      typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type  VertexOnCellIterator;

      CellContainer   cells(device.mesh());

      if (vertex_ids == 0) throw std::invalid_argument("vertices = NULL !");
      if (/*cell_id < 0 ||*/ cell_id >= cells.size()) throw std::invalid_argument("cell_id invalid !");


      VertexOnCellContainer vertices_on_cell(cells[cell_id]);
      std::size_t j = 0;
      for (VertexOnCellIterator vit = vertices_on_cell.begin(); vit != vertices_on_cell.end(); ++vit, ++j)
      {
        vertex_ids[j] = static_cast<viennashe_index_type>(vit->id().get());
      }

    } // dump_mesh


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
    typedef typename DeviceT::mesh_type              MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;

    CellContainer   cells     = viennagrid::cells(device.mesh());

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
      device.set_material(material_ids[i], cells[i]);

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
        if (doping_n[i]) device.set_doping_n(doping_n[i], cells[i]);
        if (doping_p[i]) device.set_doping_p(doping_p[i], cells[i]);
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
    typedef typename DeviceT::mesh_type              MeshType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;

    CellContainer   cells(device.mesh());

    for (std::size_t i = 0; i < len; ++i)
    {
      const viennashe_index_type id = cell_ids[i];
      if (id < cells.size())
        device.set_contact_potential(values[i], cells[id] );
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
    if (segment_id < device.segmentation().size())
      device.set_contact_potential(value, device.segment(static_cast<int>(segment_id)));
    else
      viennashe::log::warn() << "WARNING! set_contact_potential(): Invalid segment id '" << segment_id << "' ... skipping!" << std::endl;

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

    //
    // Generate device
    int_dev->stype     = libviennashe::meshtype::line_1d;
    int_dev->device_1d = new viennashe::device<viennagrid::line_1d_mesh>();

    viennashe::util::device_generation_config generator_params;
    generator_params.add_segment(0.0, len_x, static_cast<unsigned long>(points_x));
    int_dev->device_1d->generate_mesh(generator_params);

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
    if (dev != NULL)
    {
      delete (dev);
    }
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! initalize_device(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           libviennashe::initalize_device(*int_dev->device_1d,      material_ids, doping_n, doping_p); break;
    case libviennashe::meshtype::quadrilateral_2d:  libviennashe::initalize_device(*int_dev->device_quad_2d, material_ids, doping_n, doping_p); break;
    case libviennashe::meshtype::triangular_2d:     libviennashe::initalize_device(*int_dev->device_tri_2d,  material_ids, doping_n, doping_p); break;
    case libviennashe::meshtype::hexahedral_3d:     libviennashe::initalize_device(*int_dev->device_hex_3d, material_ids, doping_n, doping_p); break;
    case libviennashe::meshtype::tetrahedral_3d:    libviennashe::initalize_device(*int_dev->device_tet_3d, material_ids, doping_n, doping_p); break;
    default:
      viennashe::log::error() << "ERROR! initalize_device(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_set_material_on_segment(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           int_dev->device_1d->set_material(material_id, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::quadrilateral_2d:  int_dev->device_quad_2d->set_material(material_id, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::triangular_2d:     int_dev->device_tri_2d->set_material(material_id, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::hexahedral_3d:     int_dev->device_hex_3d->set_material(material_id, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::tetrahedral_3d:    int_dev->device_tet_3d->set_material(material_id, static_cast<int>(segment_id)); break;
    default:
      viennashe::log::error() << "ERROR! initalize_device(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_set_doping_n_on_segment(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           int_dev->device_1d->set_doping_n(doping_n, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::quadrilateral_2d:  int_dev->device_quad_2d->set_doping_n(doping_n, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::triangular_2d:     int_dev->device_tri_2d->set_doping_n(doping_n, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::hexahedral_3d:     int_dev->device_hex_3d->set_doping_n(doping_n, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::tetrahedral_3d:    int_dev->device_tet_3d->set_doping_n(doping_n, static_cast<int>(segment_id)); break;
    default:
      viennashe::log::error() << "ERROR! viennashe_set_doping_n_on_segment(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }

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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_set_doping_p_on_segment(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           int_dev->device_1d->set_doping_p(doping_p, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::quadrilateral_2d:  int_dev->device_quad_2d->set_doping_p(doping_p, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::triangular_2d:     int_dev->device_tri_2d->set_doping_p(doping_p, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::hexahedral_3d:     int_dev->device_hex_3d->set_doping_p(doping_p, static_cast<int>(segment_id)); break;
    case libviennashe::meshtype::tetrahedral_3d:    int_dev->device_tet_3d->set_doping_p(doping_p, static_cast<int>(segment_id)); break;
    default:
      viennashe::log::error() << "ERROR! viennashe_set_doping_p_on_segment(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }

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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! set_contact_potential(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           libviennashe::set_contact_potential(*int_dev->device_1d, cell_ids, values, len); break;
    case libviennashe::meshtype::quadrilateral_2d:  libviennashe::set_contact_potential(*int_dev->device_quad_2d, cell_ids, values, len); break;
    case libviennashe::meshtype::triangular_2d:     libviennashe::set_contact_potential(*int_dev->device_tri_2d, cell_ids, values, len); break;
    case libviennashe::meshtype::hexahedral_3d:     libviennashe::set_contact_potential(*int_dev->device_hex_3d, cell_ids, values, len); break;
    case libviennashe::meshtype::tetrahedral_3d:    libviennashe::set_contact_potential(*int_dev->device_tet_3d, cell_ids, values, len); break;
    default:
      viennashe::log::error() << "ERROR! initalize_device(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->device_1d == NULL)
    {
      viennashe::log::error() << "ERROR! set_contact_potential(): The device must exist!" << std::endl;
      return 1;
    }

    switch (int_dev->stype)
    {
    case libviennashe::meshtype::line_1d:           libviennashe::set_contact_potential(*int_dev->device_1d, segment_id, value); break;
    case libviennashe::meshtype::quadrilateral_2d:  libviennashe::set_contact_potential(*int_dev->device_quad_2d, segment_id, value); break;
    case libviennashe::meshtype::triangular_2d:     libviennashe::set_contact_potential(*int_dev->device_tri_2d, segment_id, value); break;
    case libviennashe::meshtype::hexahedral_3d:     libviennashe::set_contact_potential(*int_dev->device_hex_3d, segment_id, value); break;
    case libviennashe::meshtype::tetrahedral_3d:    libviennashe::set_contact_potential(*int_dev->device_tet_3d, segment_id, value); break;
    default:
      viennashe::log::error() << "ERROR! initalize_device(): UNKOWN DEVICE TYPE!" << std::endl;
      return -1;
    }
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_1d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_quad_2d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_tri_2d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_hex_3d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_tet_3d->mesh()).size());

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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_1d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_quad_2d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_tri_2d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_hex_3d->mesh()).size());
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_tet_3d->mesh()).size());
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)
      *num = static_cast<viennashe_index_type>(int_dev->device_1d->segmentation().size());
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
      *num = static_cast<viennashe_index_type>(int_dev->device_quad_2d->segmentation().size());
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)
      *num = static_cast<viennashe_index_type>(int_dev->device_tri_2d->segmentation().size());
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)
      *num = static_cast<viennashe_index_type>(int_dev->device_hex_3d->segmentation().size());
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
      *num = static_cast<viennashe_index_type>(int_dev->device_tet_3d->segmentation().size());
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


    viennashe_index_type num_seg = 0;
    viennashe_get_num_segments(dev, &num_seg);
    //
    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_1d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_quad_2d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_tri_2d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_hex_3d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::vertices(int_dev->device_tet_3d->segmentation()[static_cast<int>(segment_id)]).size());
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

    viennashe_index_type num_seg = 0;
    viennashe_get_num_segments(dev, &num_seg);
    //
    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_1d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_quad_2d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_tri_2d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_hex_3d->segmentation()[static_cast<int>(segment_id)]).size());
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)
      *num = static_cast<viennashe_index_type>(viennagrid::cells(int_dev->device_tet_3d->segmentation()[static_cast<int>(segment_id)]).size());
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)          *dim = 1;
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d) *dim = 2;
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)    *dim = 2;
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)    *dim = 3;
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)   *dim = 3;
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

    viennashe_device_impl * int_dev = (dev);

    if (int_dev->stype == libviennashe::meshtype::line_1d)          *num = 2;
    if (int_dev->stype == libviennashe::meshtype::quadrilateral_2d) *num = 4;
    if (int_dev->stype == libviennashe::meshtype::triangular_2d)    *num = 3;
    if (int_dev->stype == libviennashe::meshtype::hexahedral_3d)    *num = 6;
    if (int_dev->stype == libviennashe::meshtype::tetrahedral_3d)   *num = 4;
  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_get_num_vertices_per_cell(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device(viennashe_device * dev, viennashe_topology_type_id topology_id,
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

    viennashe::util::device_from_array_generator<viennashe_index_type> gen(vertices, cells, segmentation, num_vertices, num_cells);

    if (topology_id == viennashe_line_1d)
    {
      int_dev->stype = libviennashe::meshtype::line_1d;
      int_dev->device_1d = new viennashe::device<viennagrid::line_1d_mesh>();
      int_dev->device_1d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_quadrilateral_2d)
    {
      int_dev->stype = libviennashe::meshtype::quadrilateral_2d;
      int_dev->device_quad_2d = new viennashe::device<viennagrid::quadrilateral_2d_mesh>();
      int_dev->device_quad_2d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_triangular_2d)
    {
      int_dev->stype = libviennashe::meshtype::triangular_2d;
      int_dev->device_tri_2d = new viennashe::device<viennagrid::triangular_2d_mesh>();
      int_dev->device_tri_2d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_hexahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::hexahedral_3d;
      int_dev->device_hex_3d = new viennashe::device<viennagrid::hexahedral_3d_mesh>();
      int_dev->device_hex_3d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_tetrahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::tetrahedral_3d;
      int_dev->device_tet_3d = new viennashe::device<viennagrid::tetrahedral_3d_mesh>();
      int_dev->device_tet_3d->generate_mesh(gen);
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_create_device(): Invalid topolgy id !" << std::endl;
      return 2;
    }

    // RETURN
    *dev = int_dev;

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}



VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_flat(viennashe_device * dev, viennashe_topology_type_id topology_id,
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

    viennashe::util::device_from_flat_array_generator<viennashe_index_type> gen(vertices, cells, segmentation, num_vertices, num_cells);

    if (topology_id == viennashe_line_1d)
    {
      int_dev->stype = libviennashe::meshtype::line_1d;
      int_dev->device_1d = new viennashe::device<viennagrid::line_1d_mesh>();
      int_dev->device_1d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_quadrilateral_2d)
    {
      int_dev->stype = libviennashe::meshtype::quadrilateral_2d;
      int_dev->device_quad_2d = new viennashe::device<viennagrid::quadrilateral_2d_mesh>();
      int_dev->device_quad_2d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_triangular_2d)
    {
      int_dev->stype = libviennashe::meshtype::triangular_2d;
      int_dev->device_tri_2d = new viennashe::device<viennagrid::triangular_2d_mesh>();
      int_dev->device_tri_2d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_hexahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::hexahedral_3d;
      int_dev->device_hex_3d = new viennashe::device<viennagrid::hexahedral_3d_mesh>();
      int_dev->device_hex_3d->generate_mesh(gen);
    }
    else if (topology_id == viennashe_tetrahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::tetrahedral_3d;
      int_dev->device_tet_3d = new viennashe::device<viennagrid::tetrahedral_3d_mesh>();
      int_dev->device_tet_3d->generate_mesh(gen);
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_create_device_flat(): Invalid topolgy id !" << std::endl;
      return 2;
    }

    // RETURN
    *dev = int_dev;

  }
  catch(...)
  {
    viennashe::log::error() << "ERROR! viennashe_create_device_flat(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_create_device_from_file(viennashe_device * dev, viennashe_topology_type_id topology_id, char const * filename)
{
  try
  {
    //
    // CHECKS
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");
    CHECK_ARGUMENT_FOR_NULL(filename,2,"filename");

    viennashe_device_impl * int_dev = new viennashe_device_impl();

    if (topology_id == viennashe_line_1d)
    {
      int_dev->stype = libviennashe::meshtype::line_1d;
      int_dev->device_1d = new viennashe::device<viennagrid::line_1d_mesh>();
      int_dev->device_1d->load_mesh(std::string(filename));
    }
    else if (topology_id == viennashe_quadrilateral_2d)
    {
      int_dev->stype = libviennashe::meshtype::quadrilateral_2d;
      int_dev->device_quad_2d = new viennashe::device<viennagrid::quadrilateral_2d_mesh>();
      int_dev->device_quad_2d->load_mesh(std::string(filename));
    }
    else if (topology_id == viennashe_triangular_2d)
    {
      int_dev->stype = libviennashe::meshtype::triangular_2d;
      int_dev->device_tri_2d = new viennashe::device<viennagrid::triangular_2d_mesh>();
      int_dev->device_tri_2d->load_mesh(std::string(filename));
    }
    else if (topology_id == viennashe_hexahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::hexahedral_3d;
      int_dev->device_hex_3d = new viennashe::device<viennagrid::hexahedral_3d_mesh>();
      int_dev->device_hex_3d->load_mesh(std::string(filename));
    }
    else if (topology_id == viennashe_tetrahedral_3d)
    {
      int_dev->stype = libviennashe::meshtype::tetrahedral_3d;
      int_dev->device_tet_3d = new viennashe::device<viennagrid::tetrahedral_3d_mesh>();
      int_dev->device_tet_3d->load_mesh(std::string(filename));
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_create_device_from_file(): Invalid topolgy id !" << std::endl;
      return 2;
    }

    // RETURN
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

    // Check topology and set parameters and generate mesh
    if (dev->stype == libviennashe::meshtype::line_1d)
    {
      viennashe::util::dump_mesh(*(dev->device_1d), vertices, *num_vertices, cells, *num_cells);
    }
    else if(dev->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      viennashe::util::dump_mesh(*(dev->device_quad_2d), vertices, *num_vertices, cells, *num_cells);
    }
    else if (dev->stype == libviennashe::meshtype::triangular_2d)
    {
      viennashe::util::dump_mesh(*(dev->device_tri_2d), vertices, *num_vertices, cells, *num_cells);
    }
    else if (dev->stype == libviennashe::meshtype::hexahedral_3d)
    {
      viennashe::util::dump_mesh(*(dev->device_hex_3d), vertices, *num_vertices, cells, *num_cells);
    }
    else if (dev->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      viennashe::util::dump_mesh(*(dev->device_tet_3d), vertices, *num_vertices, cells, *num_cells);
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_grid(): Unkown topology type or malconfigured device (dev)!" << std::endl;
      return 1;
    }

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

    // Check topology and set parameters and generate mesh
    if (dev->stype == libviennashe::meshtype::line_1d)
    {
      typedef viennashe_device_impl::dev1d_type::mesh_type MeshType;
      typedef viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;

      VertexContainer vertices((dev->device_1d)->mesh());
      *x = viennagrid::point(vertices[vid])[0];
    }
    else if(dev->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      typedef viennashe_device_impl::devq2d_type::mesh_type MeshType;
      typedef viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;

      VertexContainer vertices((dev->device_quad_2d)->mesh());
      *x = viennagrid::point(vertices[vid])[0];
      *y = viennagrid::point(vertices[vid])[1];
    }
    else if (dev->stype == libviennashe::meshtype::triangular_2d)
    {
      typedef viennashe_device_impl::devt2d_type::mesh_type MeshType;
      typedef viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;

      VertexContainer vertices((dev->device_tri_2d)->mesh());
      *x = viennagrid::point(vertices[vid])[0];
      *y = viennagrid::point(vertices[vid])[1];
    }
    else if (dev->stype == libviennashe::meshtype::hexahedral_3d)
    {
      typedef viennashe_device_impl::devh3d_type::mesh_type MeshType;
      typedef viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;

      VertexContainer vertices((dev->device_hex_3d)->mesh());
      *x = viennagrid::point(vertices[vid])[0];
      *y = viennagrid::point(vertices[vid])[1];
      *z = viennagrid::point(vertices[vid])[2];
    }
    else if (dev->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      typedef viennashe_device_impl::devt3d_type::mesh_type MeshType;
      typedef viennagrid::result_of::const_vertex_range<MeshType>::type      VertexContainer;

      VertexContainer vertices((dev->device_tet_3d)->mesh());
      *x = viennagrid::point(vertices[vid])[0];
      *y = viennagrid::point(vertices[vid])[1];
      *z = viennagrid::point(vertices[vid])[2];
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): Unkown topology type or malconfigured device (dev)!" << std::endl;
      return 1;
    }

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

    // Check topology and set parameters and generate mesh
    if (dev->stype == libviennashe::meshtype::line_1d)
    {
      libviennashe::get_vertices_on_cell(*(dev->device_1d), cid, vertex_id_list);
    }
    else if(dev->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      libviennashe::get_vertices_on_cell(*(dev->device_quad_2d), cid, vertex_id_list);
    }
    else if (dev->stype == libviennashe::meshtype::triangular_2d)
    {
      libviennashe::get_vertices_on_cell(*(dev->device_tri_2d), cid, vertex_id_list);
    }
    else if (dev->stype == libviennashe::meshtype::hexahedral_3d)
    {
      libviennashe::get_vertices_on_cell(*(dev->device_hex_3d), cid, vertex_id_list);
    }
    else if (dev->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      libviennashe::get_vertices_on_cell(*(dev->device_tet_3d), cid, vertex_id_list);
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_nth_vertex(): Unkown topology type or malconfigured device (dev)!" << std::endl;
      return 1;
    }

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
