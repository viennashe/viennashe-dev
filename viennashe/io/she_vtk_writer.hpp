#ifndef VIENNASHE_IO_SHE_VTK_WRITER_HPP
#define VIENNASHE_IO_SHE_VTK_WRITER_HPP

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


// std
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

// viennagrid
#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/config.hpp"
#include "viennashe/she/postproc/all.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/io/exception.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/models/mobility.hpp"
#include "viennashe/postproc/current_density.hpp"
#include "viennashe/postproc/electric_field.hpp"
#include "viennashe/phonon/joule_heating.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/io/log_keys.h"
#include "viennashe/util/misc.hpp"


/** @file viennashe/io/she_vtk_writer.hpp
    @brief Provides a VTK writer for the computed distribution function. Outputs the (x, H)-space.
*/


namespace viennashe
{
  namespace io
  {

    /** @brief Auxiliary namespace with metafunctions */
    namespace result_of
    {
      /** @brief Meta function which translates element tags to VTK type identifiers (taking extra energy coordinate into account)
        *
        * see http://www.vtk.org/VTK/img/file-formats.pdf, Figure 2, for an overview
        *
        */
      template <typename T>
      struct she_vtk_type
      {
        //enum{ value = -1 };   //force invalid value...
        typedef typename T::ERROR_NO_SHE_VTK_WRITER_FOR_THIS_ELEMENT_TYPE_AVAILABLE  error_type;
      };

      /** @brief VTK type for a line embedded into (x, H)-space -> quadrilateral */
      template <>
      struct she_vtk_type<viennagrid::simplex_tag<1> >
      {
        enum{ value = 9 };  //VTK_quad
      };

      /** @brief VTK type for a line embedded into (x, H)-space -> quadrilateral */
      template <>
      struct she_vtk_type<viennagrid::hypercube_tag<1> >
      {
        enum{ value = 9 };  //VTK_quad
      };

      /** @brief VTK type for a quadrilateral embedded into (x, H)-space -> hexahedron */
      template <>
      struct she_vtk_type<viennagrid::quadrilateral_tag>
      {
        enum{ value = 12 };  //VTK_hexahedron
      };

      /** @brief VTK type for a triangle embedded into (x, H)-space -> wedge */
      template <>
      struct she_vtk_type<viennagrid::triangle_tag>
      {
        enum{ value = 13 };  //VTK_wedge
      };

    }


    /////////////////// VTK export ////////////////////////////

    /** @brief VTK writer class */
    template < typename SHEDeviceType,
               typename CoordSystem = typename viennagrid::result_of::coordinate_system<typename viennagrid::result_of::point<typename SHEDeviceType::mesh_type>::type>::type >
    class she_vtk_writer
    {
      protected:

      typedef typename SHEDeviceType::mesh_type                      MeshType;

      typedef typename viennagrid::result_of::cell_tag<MeshType>::type  CellTag;

      typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
      typedef typename viennagrid::result_of::vertex<MeshType>::type    VertexType;
      typedef typename viennagrid::result_of::cell<MeshType>::type      CellType;

      /** @brief Checks whether a certain cell in x-space is inside the conduction band or the valence band at total energy index index_H */
      template <typename DeviceType, typename SHEQuantity, typename CellType>
      bool is_valid(DeviceType const & device, SHEQuantity const & quan, CellType const & cell, std::size_t index_H)
      {
        if ( quan.get_unknown_index(cell, index_H) >= 0 && viennashe::materials::is_semiconductor(device.get_material(cell)) )
          return true;

        if (write_segments() && quan.get_kinetic_energy(cell, index_H) > 0 && !viennashe::materials::is_semiconductor(device.get_material(cell)))
          return true;

        return false;
      }

      /** @brief Determines the number of cells in the output mesh in (x, H)-space. */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      long get_cell_num(DeviceType const & device, SegmentType const & segment, SHEQuantity const & quan)
      {
        typedef typename viennagrid::result_of::const_vertex_range<SegmentType>::type   VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;

        typedef typename viennagrid::result_of::const_cell_range<SegmentType>::type     CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

        typedef typename viennagrid::result_of::const_vertex_range<CellType>::type      VertexOnCellContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type   VertexOnCellIterator;

        long total_cell_num = 0;

        //
        // Step 1: init write flag on vertices:
        //
        VertexContainer vertices(segment);
        vertex_write_mask_.resize(viennagrid::vertices(device.mesh()).size());
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {
          vertex_write_mask_[static_cast<std::size_t>(vit->id().get())].resize(quan.get_value_H_size());
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            vertex_write_mask_[static_cast<std::size_t>(vit->id().get())].at(index_H) = -1;
          }
        }

        //
        // Step 2: Now tag all cells where all vertices are in the conduction or valence band
        //
        CellContainer cells(segment);
        for (CellIterator cit = cells.begin();
             cit != cells.end();
            ++cit)
        {

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size()-1; ++index_H)
          {
            if (is_valid(device, quan, *cit, index_H))
            {
              ++total_cell_num;

              //tag all vertices:
              VertexOnCellContainer vertices_on_cell(*cit);
              for (VertexOnCellIterator vocit = vertices_on_cell.begin();
                  vocit != vertices_on_cell.end();
                  ++vocit)
              {
                vertex_write_mask_[static_cast<std::size_t>(vocit->id().get())].at(index_H) = 0;
                vertex_write_mask_[static_cast<std::size_t>(vocit->id().get())].at(index_H+1) = 0;
              }
            }
          }
        }

        return total_cell_num;
      }

      /** @brief Determines the number of vertices of the output mesh in (x, H)-space */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      long get_point_num(DeviceType const & device, SegmentType const & segment, SHEQuantity const & quan)
      {
        typedef typename viennagrid::result_of::const_vertex_range<SegmentType>::type   VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;

        (void)device;
        long point_num = 0;

        VertexContainer vertices(segment);
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[static_cast<std::size_t>(vit->id().get())].at(index_H) >= 0)
            {
              vertex_write_mask_[static_cast<std::size_t>(vit->id().get())].at(index_H) = point_num;
              ++point_num;
            }
          }
        }

        return point_num;
      } //writePointData



      /** @brief Writes the VTK XML file header for the unstructured grid file format */
      void writeHeader(std::ofstream & writer)
      {
        writer << "<?xml version=\"1.0\"?>" << std::endl;
        writer << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        writer << " <UnstructuredGrid>" << std::endl;
      }

      /** @brief Implementation for writing the vertex coordinates in (x, H)-space */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      void writePoints(DeviceType const & device, SegmentType const & segment, SHEQuantity const & quan, std::ofstream & writer)
      {
        typedef typename viennagrid::result_of::const_vertex_range<SegmentType>::type   VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;

        (void)device;
        writer << "   <Points>" << std::endl;
        writer << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        VertexContainer vertices(segment);
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[static_cast<std::size_t>(vit->id().get())].at(index_H) >= 0)
            {
              writer << viennagrid::point(*vit)[0] << " ";
              if (viennagrid::point(*vit).size() == 1)
                writer << "0 ";
              else
                writer << viennagrid::point(*vit)[1] << " ";
              writer << quan.get_value_H(index_H) << " ";
            }
          }
          writer << std::endl;
        }
        writer << std::endl;
        writer << "    </DataArray>" << std::endl;
        writer << "   </Points> " << std::endl;
      } //writePoints()

      /** @brief Implementation for writing the cells in (x, H)-space (derived from a mesh in x-space) */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      void writeCells(DeviceType const & device, SegmentType const & segment, SHEQuantity const & quan, std::ofstream & writer)
      {
        typedef typename viennagrid::result_of::const_cell_range<SegmentType>::type     CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

        typedef typename viennagrid::result_of::const_vertex_range<CellType>::type      VertexOnCellContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type   VertexOnCellIterator;

        writer << "   <Cells> " << std::endl;
        writer << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
        CellContainer cells(segment);

        //write prisms:
        std::size_t num_cells = 0;
        for (std::size_t index_H = 0; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          std::size_t index_H_other = index_H + 1;

          for (CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit)
          {
            if (!is_valid(device, quan, *cit, index_H))
              continue;

            VertexOnCellContainer vertices_on_cell(*cit);

            if (vertices_on_cell.size() == 2) //line segments need special treatment
            {
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[0].id().get())].at(index_H) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[1].id().get())].at(index_H) << " ";

              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[1].id().get())].at(index_H_other) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[0].id().get())].at(index_H_other) << " ";
              ++num_cells;
            }
            else if (vertices_on_cell.size() == 4) //quadrilaterals need special treatment
            {
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[0].id().get())].at(index_H) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[1].id().get())].at(index_H) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[3].id().get())].at(index_H) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[2].id().get())].at(index_H) << " ";

              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[0].id().get())].at(index_H_other) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[1].id().get())].at(index_H_other) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[3].id().get())].at(index_H_other) << " ";
              writer << vertex_write_mask_[static_cast<std::size_t>(vertices_on_cell[2].id().get())].at(index_H_other) << " ";
              ++num_cells;
            }
            else
            {
              for (VertexOnCellIterator vocit = vertices_on_cell.begin();
                  vocit != vertices_on_cell.end();
                  ++vocit)
              {
                writer << vertex_write_mask_[static_cast<std::size_t>(vocit->id().get())].at(index_H) << " ";
              }

              for (VertexOnCellIterator vocit = vertices_on_cell.begin();
                  vocit != vertices_on_cell.end();
                  ++vocit)
              {
                writer << vertex_write_mask_[static_cast<std::size_t>(vocit->id().get())].at(index_H_other) << " ";
              }

              ++num_cells;
            }
             writer << std::endl;
           }
          }

          writer << std::endl;
          writer << "    </DataArray>" << std::endl;

          writer << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
          for (std::size_t offsets = 1;
               offsets <= num_cells;
               ++offsets)
          {
            writer << (offsets * viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num * 2) << " ";
          }
          writer << std::endl;
          writer << "    </DataArray>" << std::endl;

          writer << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
          for (std::size_t offsets = 1;
                offsets <= num_cells;
                ++offsets)
          {
            writer << result_of::she_vtk_type<CellTag>::value << " ";
          }
          writer << std::endl;
          writer << "    </DataArray>" << std::endl;
          writer << "   </Cells>" << std::endl;
      }

      /** @brief Implementation for writing the data (that is usually the energy distribution function) to the vertices in (x, H)-space */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      void writePointData(DeviceType const & device,
                          SegmentType const & segment,
                          SHEQuantity const & quan,
                          std::ofstream & writer, std::string name_in_file = "result")
      {
        typedef typename viennagrid::result_of::const_vertex_range<SegmentType>::type   VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;

        (void)device;
        writer << "   <PointData Scalars=\"scalars\">" << std::endl;
        writer << "    <DataArray type=\"Float64\" Name=\"" << name_in_file << "\" format=\"ascii\">" << std::endl;

        VertexContainer vertices(segment);
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {
          for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[vit->id().get()].at(index_H) >= 0)
            {
              writer << quan.get_values(*vit, index_H)[0] << " ";
            }
          }
          writer << std::endl;
        }
        writer << std::endl;
        writer << "    </DataArray>" << std::endl;
        writer << "   </PointData>"  << std::endl;
      } //writePointData

      enum quantity_ids {
        VIENNASHE_SHE_VTK_QUAN_INVALID = 0,
        VIENNASHE_SHE_VTK_QUAN_GENERALIZED_DISTRIBUTION_FUNCTION, //generalized df
        VIENNASHE_SHE_VTK_QUAN_DISTRIBUTION_FUNCTION,
        VIENNASHE_SHE_VTK_QUAN_DENSITY_OF_STATES,
        VIENNASHE_SHE_VTK_QUAN_GROUP_VELOCITY,
        VIENNASHE_SHE_VTK_QUAN_KINETIC_ENERGY,
        VIENNASHE_SHE_VTK_QUAN_EXPANSION_ORDER,
        VIENNASHE_SHE_VTK_QUAN_EXPANSION_ADAPTION,
        VIENNASHE_SHE_VTK_QUAN_UNKNOWN_INDEX,
        VIENNASHE_SHE_VTK_QUAN_UNKNOWN_MASK,
        VIENNASHE_SHE_VTK_QUAN_UNKNOWN_NUM
      };

      /** @brief Writes data defined on cells to file */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      void writeCellDataArray(DeviceType const & device,
                              SegmentType const & segment,
                              viennashe::config const & conf,
                              SHEQuantity const & quan,
                              std::ofstream & writer,
                              quantity_ids quan_id)
      {
        typedef typename viennagrid::result_of::const_cell_range<SegmentType>::type     CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

        std::string quantity_name = quan.get_name();
        switch (quan_id)
        {
        case VIENNASHE_SHE_VTK_QUAN_GENERALIZED_DISTRIBUTION_FUNCTION: quantity_name += " (Generalized DF)"; break;
        case VIENNASHE_SHE_VTK_QUAN_DISTRIBUTION_FUNCTION:             quantity_name += " (DF)"; break;
        case VIENNASHE_SHE_VTK_QUAN_DENSITY_OF_STATES:                 quantity_name += " (density of states)"; break;
        case VIENNASHE_SHE_VTK_QUAN_GROUP_VELOCITY:                    quantity_name += " (group velocity)"; break;
        case VIENNASHE_SHE_VTK_QUAN_KINETIC_ENERGY:                    quantity_name += " (kinetic energy)"; break;
        case VIENNASHE_SHE_VTK_QUAN_EXPANSION_ORDER:                   quantity_name += " (expansion order)"; break;
        case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_INDEX:                     quantity_name += " (unknown index)"; break;
        case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_MASK:                      quantity_name += " (unknown mask)"; break;
        case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_NUM:                       quantity_name += " (unknown number)"; break;
        default: throw std::runtime_error("Internal error: Unknown quan_id in she_vtk_writer::writeCellDataArray()");
        }

        typename viennashe::config::dispersion_relation_type dispersion = conf.dispersion_relation(quan.get_carrier_type_id());

        writer << "    <DataArray type=\"Float64\" Name=\"Generalized " << quantity_name << "\" format=\"ascii\">" << std::endl;

        CellContainer cells(segment);

        //write prisms:
        for (std::size_t index_H = 0; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          double value = 0;
          double dos;

          for (CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit)
          {
            if (!is_valid(device, quan, *cit, index_H))
              continue;

            switch (quan_id)
            {
            case VIENNASHE_SHE_VTK_QUAN_GENERALIZED_DISTRIBUTION_FUNCTION:
              switch (conf.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                dos = averaged_density_of_states(quan, dispersion, *cit, index_H);
                value = quan.get_values(*cit, index_H)[0] * dos;
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                value = quan.get_values(*cit, index_H)[0];
                break;
              }
              break;

            case VIENNASHE_SHE_VTK_QUAN_DISTRIBUTION_FUNCTION:
              switch (conf.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                value = quan.get_values(*cit, index_H)[0];
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                dos = averaged_density_of_states(quan, dispersion, *cit, index_H);
                value = (dos > 0) ? quan.get_values(*cit, index_H)[0] / dos : 0;  //constant continuation of DF towards zero
                break;
              }
              break;
            case VIENNASHE_SHE_VTK_QUAN_DENSITY_OF_STATES:                 value = dispersion.density_of_states(quan.get_kinetic_energy(*cit, index_H)); break;
            case VIENNASHE_SHE_VTK_QUAN_GROUP_VELOCITY:                    value = dispersion.velocity(quan.get_kinetic_energy(*cit, index_H)); break;
            case VIENNASHE_SHE_VTK_QUAN_KINETIC_ENERGY:                    value = quan.get_kinetic_energy(*cit, index_H); break;
            case VIENNASHE_SHE_VTK_QUAN_EXPANSION_ORDER:                   value = static_cast<double>(quan.get_expansion_order(*cit, index_H)); break;
            case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_INDEX:                     value = static_cast<double>(quan.get_unknown_index(*cit, index_H)); break;
            case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_MASK:                      value = static_cast<double>(quan.get_unknown_mask(*cit, index_H)); break;
            case VIENNASHE_SHE_VTK_QUAN_UNKNOWN_NUM:                       value = static_cast<double>(quan.get_unknown_num(*cit, index_H)); break;
            default: throw std::runtime_error("Internal error: Unknown quan_id in she_vtk_writer::writeCellDataArray()");
            }

            writer << value << " ";
          }
          writer << std::endl;
        }
        writer << "    </DataArray>" << std::endl;
      } //writeCellDataArray

      /** @brief Writes data defined on cells to file */
      template <typename DeviceType, typename SegmentType, typename SHEQuantity>
      void writeCellData(DeviceType const & device,
                         SegmentType const & segment,
                         viennashe::config const & conf,
                         SHEQuantity const & quan,
                         std::ofstream & writer)
      {
        writer << "   <CellData Scalars=\"scalars\">" << std::endl;

        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_GENERALIZED_DISTRIBUTION_FUNCTION);
        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_DISTRIBUTION_FUNCTION);
        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_DENSITY_OF_STATES);
        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_GROUP_VELOCITY);
        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_KINETIC_ENERGY);
        writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_EXPANSION_ORDER);
        if (with_debug_quantities())
        {
          writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_INDEX);
          writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_MASK);
          writeCellDataArray(device, segment, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_NUM);
        }

        writer << "   </CellData>"  << std::endl;
      } //writeCellData

      void writeFooter(std::ofstream & writer)
      {
        writer << " </UnstructuredGrid>" << std::endl;
        writer << "</VTKFile>" << std::endl;
      }

      template <typename DeviceType, typename SegmentType>
      bool segment_is_semiconductor_only(DeviceType const & device, SegmentType const & segment)
      {
        typedef typename viennagrid::result_of::const_cell_range<SegmentType>::type    CellSegmentContainer;
        typedef typename viennagrid::result_of::iterator<CellSegmentContainer>::type   CellSegmentIterator;

        CellSegmentContainer cells(segment);
        for (CellSegmentIterator cit  = cells.begin();
                                 cit != cells.end();
                               ++cit)
        {
          if (!viennashe::materials::is_semiconductor(device.get_material(*cit)))
            return false;
        }
        return true;
      }

      template <typename DeviceType, typename SegmentType, typename SHEQuantityT>
      void write_segment(DeviceType const & device,
                         SegmentType const & segment,
                         viennashe::config const & conf,
                         SHEQuantityT const & quan,
                         std::string filename)
      {
        std::ofstream writer(filename.c_str());

        if (!writer)
          throw cannot_open_file_exception(filename);

        writeHeader(writer);

        long cell_num = get_cell_num(device, segment, quan); //important: get_cell_num() prior to get_point_num()!!
        long point_num = get_point_num(device, segment, quan);

        writer << "  <Piece NumberOfPoints=\""
              << point_num
              << "\" NumberOfCells=\""
              << cell_num
              << "\">" << std::endl;


        writePoints(device, segment, quan, writer);
        if ( (write_segments_ && segment_is_semiconductor_only(device, segment)) || !write_segments_)
          writeCellData(device, segment, conf, quan, writer);
        writeCells(device, segment, quan, writer);

        writer << "  </Piece>" << std::endl;

        writeFooter(writer);

      }

    public:

      she_vtk_writer() : write_segments_(false), with_debug_quantities_(false) {}

      /** @brief Triggers the write process
       *
       * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
       * @param quan           The SHE quantity in (x, H)-space to be written (typically distribution function, expansion order or error indicator)
       * @param filename       Name of the file to be written to
       * @param conf           The simulator configuration
       */
      template <typename DeviceType, typename SHEQuantityT>
      void operator()(DeviceType const & device,
                      viennashe::config const & conf,
                      SHEQuantityT const & quan,
                      std::string const & filename)
      {
        typedef typename viennagrid::result_of::segmentation<MeshType>::type    SegmentationType;
        typedef typename SegmentationType::segment_handle_type                  SegmentHandleType;

        if (write_segments_)
        {
          //
          // Step 1: Write meta information (pvd file)
          //
          std::stringstream ss;
          ss << filename << "_main.pvd";
          std::ofstream writer(ss.str().c_str());

          std::string short_filename = filename;
          std::string::size_type pos = filename.rfind("/");
          if (pos == std::string::npos)
            pos = filename.rfind("\\");   //A tribute to Windows

          if (pos != std::string::npos)
            short_filename = filename.substr(pos+1, filename.size());

          writer << "<?xml version=\"1.0\"?>" << std::endl;
          writer << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
          writer << "<Collection>" << std::endl;

          for (typename SegmentationType::const_iterator it = device.segmentation().begin(); it != device.segmentation().end(); ++it)
          {
            SegmentHandleType const & seg = *it;
            writer << "    <DataSet part=\"" << seg.id() << "\" file=\"" << short_filename << "_" << seg.id() << ".vtu\" name=\"Segment_" << seg.id() << "\"/>" << std::endl;
          }

          writer << "  </Collection>" << std::endl;
          writer << "</VTKFile>" << std::endl;
          writer.close();

          //
          // Step 2: Write each segment to a separate .vtu file:
          //
          for (typename SegmentationType::const_iterator it = device.segmentation().begin(); it != device.segmentation().end(); ++it)
          {
            SegmentHandleType const & seg = *it;

            ss.str("");
            ss << filename << "_" << seg.id() << ".vtu";
            write_segment(device, *it, conf, quan, ss.str());
          }
        }
        else
        {
          //
          // Write full mesh to a single .vtu file:
          //
          write_segment(device, device.mesh(), conf, quan, filename + ".vtu");
        }

      }

      bool write_segments() const { return write_segments_; }
      void write_segments(bool b) { write_segments_ = b; }

      bool with_debug_quantities() const { return with_debug_quantities_; }
      void with_debug_quantities(bool b) { with_debug_quantities_ = b; }

    private:
      std::vector<viennashe::she_index_vector_type>  vertex_write_mask_;
      bool write_segments_;
      bool with_debug_quantities_;

    }; //she_vtk_writer


    ///////////////////// Convenience routines /////////////////////////////


    /** @brief Compatibility overload in order to enable compilation for three-dimensional devices. No action if called. */
    template < typename DeviceType >
    class she_vtk_writer<DeviceType, viennagrid::cartesian_cs<3> >
    {
      public:
        template <typename DataType>
        void operator()(DeviceType const & device,
                        DataType const & data,
                        std::string const & filename, std::string name_in_file = "result")
        {
          (void)device; (void)data; (void)filename; (void)name_in_file;
          log::warning() << "* she_vtk_writer::operator(): Cannot write distribution function for 3d devices. Skipping..." << std::endl;
        }
    };

    /** @brief Convenience routine for writing a single macroscopic quantity to a VTK file.
     *
     * @param quantity     An accessor for a macroscopic quantity
     * @param device       The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param filename     Name of the file to be written to
     * @param name_in_file   The quantity name to be used in the VTK file
     */
    template <typename QuantityType,
              typename DeviceType>
    void write_vertex_quantity_to_VTK_file(QuantityType const & quantity,
                                           DeviceType const & device,
                                           std::string filename,
                                           std::string name_in_file = "viennashe_quantity")
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::vertex<MeshType>::type                 VertexType;
      typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type     VertexContainer;
      typedef typename viennagrid::result_of::iterator<VertexContainer>::type        VertexIterator;

      MeshType const & mesh = device.mesh();

      VertexContainer vertices(mesh);
      std::vector<double> vtk_data(vertices.size());
      for (VertexIterator vit = vertices.begin();
          vit != vertices.end();
          ++vit)
      {
        vtk_data[vit->id().get()] = quantity(*vit);
      }

      log::info<log_she_vtk_writer>() << "* write_quantity_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. ParaView)" << std::endl;

      viennagrid::io::vtk_writer<MeshType> my_vtk_writer;
      my_vtk_writer.add_scalar_data_on_vertices(viennagrid::make_accessor<VertexType>(vtk_data), name_in_file);
      my_vtk_writer(mesh, device.segmentation(), filename);
    }

    template <typename DeviceType>
    void write_cell_quantity_to_VTK_file(std::vector<double> const & vtk_data,
                                         DeviceType const & device,
                                         std::string filename,
                                         std::string name_in_file = "viennashe_quantity")
    {
      typedef typename DeviceType::mesh_type                                       MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                 CellType;

      log::info<log_she_vtk_writer>() << "* write_quantity_to_VTK_file(): Writing data to '" << filename
                                      << "' (can be viewed with e.g. ParaView)" << std::endl;

      viennagrid::io::vtk_writer<MeshType> my_vtk_writer;
      my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(vtk_data), name_in_file);
      my_vtk_writer(device.mesh(), device.segmentation(), filename);
    }

    /** @brief Convenience routine for writing a single macroscopic quantity to a VTK file.
     *
     * @param quantity     An accessor for a macroscopic quantity
     * @param device       The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param filename     Name of the file to be written to
     * @param name_in_file   The quantity name to be used in the VTK file
     */
    template <typename QuantityType,
              typename DeviceType>
    void write_cell_quantity_to_VTK_file(QuantityType const & quantity,
                                         DeviceType const & device,
                                         std::string filename,
                                         std::string name_in_file = "viennashe_quantity")
    {
      typedef typename DeviceType::mesh_type                                       MeshType;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type        CellIterator;

      CellContainer cells(device.mesh());
      std::vector<double> vtk_data(cells.size());
      for (CellIterator cit = cells.begin();
           cit != cells.end();
           ++cit)
      {
        vtk_data[static_cast<std::size_t>(cit->id().get())] = quantity(*cit);
      }

      write_cell_quantity_to_VTK_file(vtk_data, device, filename, name_in_file);
    }

    namespace result_of
    {
      /** @brief Helper routine for extracting the ViennaGrid topology tag from a quantity wrapper. Works well with viennashe::util::spatial_quantity_wrapper. Custom wrappers should add template specializations of this helper metafunction */
      template <typename T>
      struct topology_tag //by default, every quantity is assumed to write to a vertex.
      {
        typedef viennagrid::vertex_tag type;
      };

      template <typename DeviceType, typename TagT>
      struct topology_tag< viennashe::util::spatial_quantity_wrapper<DeviceType, TagT> >
      {
        typedef TagT type;
      };

    }

    /** @brief Namespace for implementation details within viennashe::io. Typically not of interest for a library user. */
    namespace detail
    {
      template <typename QuantityType,
                typename DeviceType>
      void write_quantity_to_VTK_file(QuantityType const & quantity,
                                      DeviceType const & device,
                                      std::string filename,
                                      std::string name_in_file,
                                      viennagrid::vertex_tag)
      {
        write_cell_quantity_to_VTK_file(quantity, device, filename, name_in_file);
      }

      template <typename QuantityType,
                typename DeviceType,
                typename Tag>
      void write_quantity_to_VTK_file(QuantityType const & quantity,
                                      DeviceType const & device,
                                      std::string filename,
                                      std::string name_in_file,
                                      Tag)
      {
        write_cell_quantity_to_VTK_file(quantity, device, filename, name_in_file);
      }
    }

    /** @brief Generic interface function for writing a quantity to a VTK file. Automatically dispatches between vertex and cell quantity.
     *
     *    Custom functors may need to overload detail::extract_topology_tag if writing cell quantities (by default, unidentified quantities are assumed to be vertex-quantities)
     */
    template <typename QuantityType,
              typename DeviceType>
    void write_quantity_to_VTK_file(QuantityType const & quantity,
                                    DeviceType const & device,
                                    std::string filename,
                                    std::string name_in_file = "viennashe_quantity")
    {
      write_cell_quantity_to_VTK_file(quantity, device, filename, name_in_file);
    }



    namespace detail
    {
      template <typename ContainerType, typename KeyType, typename ValueType>
      class container_accessor
      {
      public:
        typedef ValueType  value_type;

        container_accessor(ContainerType const & container) : container_(container) {}

        value_type const & operator()(KeyType const & key) const { cached_value_ = container_.at(static_cast<std::size_t>(key.id().get())); return cached_value_; }
        value_type const & operator[](KeyType const & key) const { cached_value_ = container_.at(static_cast<std::size_t>(key.id().get())); return cached_value_; }
        value_type const &         at(KeyType const & key) const { cached_value_ = container_.at(static_cast<std::size_t>(key.id().get())); return cached_value_; }

        value_type const * find(KeyType const &) const { return NULL; }
        /*{
          return &(container_[key.id().get()]);
        }*/

      private:
        ContainerType const & container_;
        mutable value_type cached_value_;
      };

      template <typename KeyType, typename ValueType, typename ContainerType>
      container_accessor<ContainerType, KeyType, ValueType> make_accessor(ContainerType const & c)
      {
        return container_accessor<ContainerType, KeyType, ValueType>(c);
      }

    } // namespace detail

    /** @brief Generic interface function for writing simulated quantities to a VTK file. */
    template <typename DeviceType>
    void write_quantities_to_VTK_file(viennashe::simulator<DeviceType> const & simulator_obj,
                                      std::string filename,
                                      bool include_debug_information = false)
    {
      typedef typename viennashe::simulator<DeviceType>               SimulatorType;
      typedef typename DeviceType::mesh_type                          MeshType;
      typedef typename viennagrid::result_of::cell<MeshType>::type    CellType;

      typedef typename SimulatorType::unknown_quantity_type    UnknownQuantityType;

      DeviceType const & device = simulator_obj.device();

      viennagrid::io::vtk_writer<MeshType> my_vtk_writer;

      std::deque<std::vector<double> > electric_field(viennagrid::cells(device.mesh()).size(), std::vector<double>(3));
      std::deque<std::vector<double> > current_n(viennagrid::cells(device.mesh()).size(), std::vector<double>(3));
      std::deque<std::vector<double> > current_p(viennagrid::cells(device.mesh()).size(), std::vector<double>(3));
      std::deque<std::vector<double> > carrier_velocity_n(viennagrid::cells(device.mesh()).size(), std::vector<double>(3));
      std::deque<std::vector<double> > carrier_velocity_p(viennagrid::cells(device.mesh()).size(), std::vector<double>(3));
      std::deque<double> pwr_density_container(viennagrid::cells(device.mesh()).size());
      std::deque<double> avg_energy_n(viennagrid::cells(device.mesh()).size());
      std::deque<double> avg_energy_p(viennagrid::cells(device.mesh()).size());
      std::deque<double> avg_trap_occupancy(viennagrid::cells(device.mesh()).size());

      //
      // Device data
      //
      my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(device.doping_n()), "Donator Doping");
      my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(device.doping_p()), "Acceptor Doping");
      my_vtk_writer.add_scalar_data_on_cells(    detail::make_accessor<CellType, double>(device.material()), "Material IDs");


      //
      // Simulator data
      //
      std::deque<UnknownQuantityType> const & unknown_quans = simulator_obj.quantity_history(0).unknown_quantities();

      for (std::size_t quan_index = 0; quan_index < unknown_quans.size(); ++quan_index)
      {
        UnknownQuantityType const & quan = unknown_quans.at(quan_index);

        bool quantity_is_from_she = false;

        // electric field
        if (quan.get_name() == viennashe::quantity::potential())
        {
          viennashe::write_electric_field_to_container(device,
                                                       simulator_obj.quantities().potential(),
                                                       electric_field);
          my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(electric_field), "Electric Field");
        }

        // electron current
        if (quan.get_name() == viennashe::quantity::electron_density())
        {
          if (simulator_obj.config().get_electron_equation() == viennashe::EQUATION_SHE)
          {
            quantity_is_from_she = true;
            viennashe::she::write_current_density_to_container(device,
                                                               simulator_obj.config(),
                                                               simulator_obj.quantities().electron_distribution_function(),
                                                               current_n);
            viennashe::she::write_carrier_velocity_to_container(device,
                                                               simulator_obj.config(),
                                                               simulator_obj.quantities().electron_distribution_function(),
                                                               carrier_velocity_n);
            my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(carrier_velocity_n), "Electron avg. carrier drift velocity");
            viennashe::she::write_kinetic_carrier_energy_to_container(device,
                                                                      simulator_obj.config(),
                                                                      simulator_obj.quantities().electron_distribution_function(),
                                                                      avg_energy_n);
            my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(avg_energy_n), "Electron avg. energy");
          }
          else
          {
            viennashe::write_current_density_to_container(simulator_obj,
                                                          viennashe::ELECTRON_TYPE_ID,
                                                          current_n);
          }
          my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(current_n), "Electron current density");
        }

        // hole current
        if (quan.get_name() == viennashe::quantity::hole_density())
        {
          if (simulator_obj.config().get_hole_equation() == viennashe::EQUATION_SHE)
          {
            quantity_is_from_she = true;
            viennashe::she::write_current_density_to_container(device,
                                                               simulator_obj.config(),
                                                               simulator_obj.quantities().hole_distribution_function(),
                                                               current_p);
            viennashe::she::write_carrier_velocity_to_container(device,
                                                               simulator_obj.config(),
                                                               simulator_obj.quantities().hole_distribution_function(),
                                                               carrier_velocity_p);
            my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(carrier_velocity_p), "Hole avg. carrier drift velocity");
            viennashe::she::write_kinetic_carrier_energy_to_container(device,
                                                                      simulator_obj.config(),
                                                                      simulator_obj.quantities().hole_distribution_function(),
                                                                      avg_energy_p);
            my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(avg_energy_p), "Hole avg. energy");
          }
          else
          {
            viennashe::write_current_density_to_container(simulator_obj,
                                                          viennashe::HOLE_TYPE_ID,
                                                          current_p);
          }
          my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(current_p), "Hole current density");
        }


        // lattice temperature
        /*
        if (quan.get_name() == viennashe::quantity::lattice_temperature())
        {
          typedef typename SimulatorType::SHETimeStepQuantitiesT QuantitiesType;
          typedef typename viennashe::hde::power_density_accessor<DeviceType, QuantitiesType> PowerDensityAccessorType;

          PowerDensityAccessorType pdacc(device,
                                         simulator_obj.quantities(),
                                         simulator_obj.config());

          viennashe::write_macroscopic_quantity_to_container(device, pdacc, pwr_density_container);

          my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(pwr_density_container), "Joule heating power density");
        } */


        my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(quan.values()), quan.get_name());
        if (!quantity_is_from_she && include_debug_information)
        {
          my_vtk_writer.add_scalar_data_on_cells(detail::make_accessor<CellType, double>(quan.boundary_types()), quan.get_name() + " boundary types");
          my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(quan.boundary_values()), quan.get_name() + " boundary values");
          my_vtk_writer.add_scalar_data_on_cells(detail::make_accessor<CellType, double>(quan.defined_but_unknown_mask()), quan.get_name() + " mask");
          my_vtk_writer.add_scalar_data_on_cells(detail::make_accessor<CellType, double>(quan.unknowns_indices()), quan.get_name() + " indices");
        }
      } // for unkown quans


      if (simulator_obj.config().with_traps())
      {
        typedef typename viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type        CellIterator;

        CellContainer cells(device.mesh());
        for (CellIterator cit  = cells.begin();
                          cit != cells.end();
                        ++cit)
        {
          avg_trap_occupancy[static_cast<std::size_t>(cit->id().get())] = 0;
          // Average ...
          double ft = 0.0;

          const std::size_t num_traps =  simulator_obj.quantities().num_trap_unknown_indices(*cit);

          if (num_traps == 0) continue;

          for (std::size_t i = 0; i < num_traps; ++i)
          {
            ft += simulator_obj.quantities().trap_occupancy(*cit, i);
          }

          ft = ft / num_traps;
          avg_trap_occupancy[static_cast<std::size_t>(cit->id().get())] = ft;

        }
        my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(avg_trap_occupancy), "Average trap occupancy");
      } // traps

      my_vtk_writer(device.mesh(), device.segmentation(), filename);

    }

  } //namespace io
} //namespace viennashe

#endif
