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
#include "viennagrid/viennagrid.h"

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

    inline int viennagrid_to_vtk_type(viennagrid_element_type type)
    {
      switch (type)
      {
      case VIENNAGRID_ELEMENT_TYPE_LINE:          return  9; // quadrilateral
      case VIENNAGRID_ELEMENT_TYPE_TRIANGLE:      return 13; // wedge
      case VIENNAGRID_ELEMENT_TYPE_QUADRILATERAL: return 12; // hexahedron
      default:
        throw std::runtime_error("ViennaGrid element type not recognized.");
      }
    }


    /////////////////// VTK export ////////////////////////////

    /** @brief VTK writer class */
    template<typename SHEDeviceType>
    class she_vtk_writer
    {
    protected:

      typedef typename SHEDeviceType::mesh_type                      MeshType;

      /** @brief Checks whether a certain cell in x-space is inside the conduction band or the valence band at total energy index index_H */
      template <typename DeviceType, typename SHEQuantity>
      bool is_valid(DeviceType const & device, SHEQuantity const & quan, viennagrid_element_type cell, std::size_t index_H)
      {
        if ( quan.get_unknown_index(cell, index_H) >= 0 && viennashe::materials::is_semiconductor(device.get_material(cell)) )
          return true;

        if (write_segments() && quan.get_kinetic_energy(cell, index_H) > 0 && !viennashe::materials::is_semiconductor(device.get_material(cell)))
          return true;

        return false;
      }

      /** @brief Determines the number of cells in the output mesh in (x, H)-space. */
      template <typename DeviceType, typename SHEQuantity>
      long get_cell_num(DeviceType const & device, SHEQuantity const & quan)
      {

        long total_cell_num = 0;

        //
        // Step 1: init write flag on vertices:
        //
        viennagrid_element_id *vertices_begin, *vertices_end;
        viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

        vertex_write_mask_.resize(vertices_end - vertices_begin);

        for (viennagrid_element_id *vit  = vertices_begin;
                                    vit != vertices_end;
                                  ++vit)
        {
          std::size_t index(viennagrid_index_from_element_id(*vit));

          vertex_write_mask_[index].resize(quan.get_value_H_size());
          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            vertex_write_mask_[index].at(index_H) = -1;
          }
        }

        //
        // Step 2: Now tag all cells where all vertices are in the conduction or valence band
        //
        viennagrid_dimension cell_dim;
        viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

        viennagrid_element_id *cells_begin, *cells_end;
        viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size()-1; ++index_H)
          {
            if (is_valid(device, quan, *cit, index_H))
            {
              ++total_cell_num;

              //tag all vertices:
              viennagrid_element_id *vertices_on_cell_begin, *vertices_on_cell_end;
              viennagrid_element_boundary_elements(device.mesh(), *cit, 0, &vertices_on_cell_begin, &vertices_on_cell_end);

              for (viennagrid_element_id *vocit  = vertices_on_cell_begin;
                                          vocit != vertices_on_cell_end;
                                        ++vocit)
              {
                std::size_t index(viennagrid_index_from_element_id(*vocit));

                vertex_write_mask_[index].at(index_H) = 0;
                vertex_write_mask_[index].at(index_H+1) = 0;
              }
            }
          }
        }

        return total_cell_num;
      }

      /** @brief Determines the number of vertices of the output mesh in (x, H)-space */
      template <typename DeviceType, typename SHEQuantity>
      long get_point_num(DeviceType const & device, SHEQuantity const & quan)
      {
        long point_num = 0;

        viennagrid_element_id *vertices_begin, *vertices_end;
        viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

        for (viennagrid_element_id *vit  = vertices_begin;
                                    vit != vertices_end;
                                  ++vit)
        {
          std::size_t index(viennagrid_index_from_element_id(*vit));

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[index].at(index_H) >= 0)
            {
              vertex_write_mask_[index].at(index_H) = point_num;
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
      template <typename DeviceType, typename SHEQuantity>
      void writePoints(DeviceType const & device, SHEQuantity const & quan, std::ofstream & writer)
      {
        writer << "   <Points>" << std::endl;
        writer << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        viennagrid_element_id *vertices_begin, *vertices_end;
        viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

        for (viennagrid_element_id *vit  = vertices_begin;
                                    vit != vertices_end;
                                  ++vit)
        {
          viennagrid_numeric *coords;
          viennagrid_mesh_vertex_coords_get(device.mesh(), *vit, &coords);

          viennagrid_dimension geo_dim;
          viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim);

          for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(*vit))].at(index_H) >= 0)
            {
              writer << coords[0] << " ";
              if (geo_dim == 1)
                writer << "0 ";
              else
                writer << coords[1] << " ";
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
      template <typename DeviceType, typename SHEQuantity>
      void writeCells(DeviceType const & device, SHEQuantity const & quan, std::ofstream & writer)
      {
        writer << "   <Cells> " << std::endl;
        writer << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;

        viennagrid_dimension cell_dim;
        viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

        viennagrid_element_id *cells_begin, *cells_end;
        viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

        std::vector<std::size_t> offsets(1);
        std::vector<std::size_t> types;

        //write prisms:
        std::size_t num_cells = 0;
        for (std::size_t index_H = 0; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          std::size_t index_H_other = index_H + 1;

          for (viennagrid_element_id *cit  = cells_begin;
                                      cit != cells_end;
                                    ++cit)
          {
            if (!is_valid(device, quan, *cit, index_H))
              continue;

            viennagrid_element_id *vertices_on_cell_begin, *vertices_on_cell_end;
            viennagrid_element_boundary_elements(device.mesh(), *cit, 0, &vertices_on_cell_begin, &vertices_on_cell_end);

            offsets.push_back(offsets.back() + 2 * (vertices_on_cell_end - vertices_on_cell_begin));

            viennagrid_element_type element_type;
            viennagrid_element_type_get(device.mesh(), *cit, &element_type);
            types.push_back(viennagrid_to_vtk_type(element_type));

            if (vertices_on_cell_end - vertices_on_cell_begin == 2) //line segments need special treatment
            {
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[0]))].at(index_H) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[1]))].at(index_H) << " ";

              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[1]))].at(index_H_other) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[0]))].at(index_H_other) << " ";
              ++num_cells;
            }
            else if (vertices_on_cell_end - vertices_on_cell_begin == 4) //quadrilaterals need special treatment
            {
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[0]))].at(index_H) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[1]))].at(index_H) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[3]))].at(index_H) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[2]))].at(index_H) << " ";

              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[0]))].at(index_H_other) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[1]))].at(index_H_other) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[3]))].at(index_H_other) << " ";
              writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(vertices_on_cell_begin[2]))].at(index_H_other) << " ";
              ++num_cells;
            }
            else
            {
              for (viennagrid_element_id *vocit  = vertices_on_cell_begin;
                                          vocit != vertices_on_cell_end;
                                        ++vocit)
              {
                writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(*vocit))].at(index_H) << " ";
              }

              for (viennagrid_element_id *vocit  = vertices_on_cell_begin;
                                          vocit != vertices_on_cell_end;
                                        ++vocit)
              {
                writer << vertex_write_mask_[std::size_t(viennagrid_index_from_element_id(*vocit))].at(index_H_other) << " ";
              }

              ++num_cells;
            }
             writer << std::endl;
           }
          }

          writer << std::endl;
          writer << "    </DataArray>" << std::endl;

          writer << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
          for (std::size_t i = 1; i < offsets.size(); ++i)
            writer << offsets[i] << " ";

          writer << std::endl;
          writer << "    </DataArray>" << std::endl;

          writer << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
          for (std::size_t i = 0; i < types.size(); ++i)
            writer << types[i] << " ";

          writer << std::endl;
          writer << "    </DataArray>" << std::endl;
          writer << "   </Cells>" << std::endl;
      }

      /** @brief Implementation for writing the data (that is usually the energy distribution function) to the vertices in (x, H)-space */
      template <typename DeviceType, typename SHEQuantity>
      void writePointData(DeviceType const & device,
                          SHEQuantity const & quan,
                          std::ofstream & writer, std::string name_in_file = "result")
      {
        writer << "   <PointData Scalars=\"scalars\">" << std::endl;
        writer << "    <DataArray type=\"Float64\" Name=\"" << name_in_file << "\" format=\"ascii\">" << std::endl;

        viennagrid_element_id *vertices_begin, *vertices_end;
        viennagrid_mesh_elements_get(device.mesh(), 0, &vertices_begin, &vertices_end);

        for (viennagrid_element_id *vit  = vertices_begin;
                                    vit != vertices_end;
                                  ++vit)
        {
          for (size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
          {
            if (vertex_write_mask_[viennagrid_index_from_element_id(*vit)].at(index_H) >= 0)
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
      template <typename DeviceType, typename SHEQuantity>
      void writeCellDataArray(DeviceType const & device,
                              viennashe::config const & conf,
                              SHEQuantity const & quan,
                              std::ofstream & writer,
                              quantity_ids quan_id)
      {
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

        viennagrid_dimension cell_dim;
        viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

        viennagrid_element_id *cells_begin, *cells_end;
        viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);

        //write prisms:
        for (std::size_t index_H = 0; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          double value = 0;
          double dos;

          for (viennagrid_element_id *cit  = cells_begin;
                                      cit != cells_end;
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
      template <typename DeviceType, typename SHEQuantity>
      void writeCellData(DeviceType const & device,
                         viennashe::config const & conf,
                         SHEQuantity const & quan,
                         std::ofstream & writer)
      {
        writer << "   <CellData Scalars=\"scalars\">" << std::endl;

        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_GENERALIZED_DISTRIBUTION_FUNCTION);
        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_DISTRIBUTION_FUNCTION);
        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_DENSITY_OF_STATES);
        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_GROUP_VELOCITY);
        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_KINETIC_ENERGY);
        writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_EXPANSION_ORDER);
        if (with_debug_quantities())
        {
          writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_INDEX);
          writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_MASK);
          writeCellDataArray(device, conf, quan, writer, VIENNASHE_SHE_VTK_QUAN_UNKNOWN_NUM);
        }

        writer << "   </CellData>"  << std::endl;
      } //writeCellData

      void writeFooter(std::ofstream & writer)
      {
        writer << " </UnstructuredGrid>" << std::endl;
        writer << "</VTKFile>" << std::endl;
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

        long cell_num = get_cell_num(device, quan); //important: get_cell_num() prior to get_point_num()!!
        long point_num = get_point_num(device, quan);

        writer << "  <Piece NumberOfPoints=\""
              << point_num
              << "\" NumberOfCells=\""
              << cell_num
              << "\">" << std::endl;


        writePoints(device, quan, writer);
        //if ( (write_segments_ && segment_is_semiconductor_only(device, segment)) || !write_segments_)
          writeCellData(device, conf, quan, writer);
        writeCells(device, quan, writer);

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

        /*if (write_segments_)
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
        {*/
          //
          // Write full mesh to a single .vtu file:
          //
          write_segment(device, device.mesh(), conf, quan, filename + ".vtu");
        //}

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



    /** @brief Convenience routine for writing a single macroscopic quantity to a VTK file.
     *
     * @param quantity     An accessor for a macroscopic quantity
     * @param device       The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param filename     Name of the file to be written to
     * @param name_in_file   The quantity name to be used in the VTK file
     */
    /* TODO: Migrate to ViennaGrid 3.0
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
    } */


    /** @brief Convenience routine for writing a single macroscopic quantity to a VTK file.
     *
     * @param quantity     An accessor for a macroscopic quantity
     * @param device       The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param filename     Name of the file to be written to
     * @param name_in_file   The quantity name to be used in the VTK file
     */
    /* TODO: Migrate to ViennaGrid 3.0
    template <typename QuantityType,
              typename DeviceType>
    void write_cell_quantity_to_VTK_file(QuantityType const & quantity,
                                         DeviceType const & device,
                                         std::string filename,
                                         std::string name_in_file = "viennashe_quantity")
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                 CellType;
      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type     CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type        CellIterator;

      MeshType const & mesh = device.mesh();

      CellContainer cells(mesh);
      std::vector<double> vtk_data(cells.size());
      for (CellIterator cit = cells.begin();
           cit != cells.end();
           ++cit)
      {
        vtk_data[static_cast<std::size_t>(cit->id().get())] = quantity(*cit);
      }

      log::info<log_she_vtk_writer>() << "* write_quantity_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. ParaView)" << std::endl;

      viennagrid::io::vtk_writer<MeshType> my_vtk_writer;
      my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(vtk_data), name_in_file);
      my_vtk_writer(mesh, device.segmentation(), filename);
    }*/


    /** @brief Namespace for implementation details within viennashe::io. Typically not of interest for a library user. */
    /*namespace detail
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
    }*/

    /** @brief Generic interface function for writing a quantity to a VTK file. Automatically dispatches between vertex and cell quantity.
     *
     *    Custom functors may need to overload detail::extract_topology_tag if writing cell quantities (by default, unidentified quantities are assumed to be vertex-quantities)
     */
    /*
    template <typename QuantityType,
              typename DeviceType>
    void write_quantity_to_VTK_file(QuantityType const & quantity,
                                    DeviceType const & device,
                                    std::string filename,
                                    std::string name_in_file = "viennashe_quantity")
    {
      write_cell_quantity_to_VTK_file(quantity, device, filename, name_in_file);
    } */


/*
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
*/
    /** @brief Generic interface function for writing simulated quantities to a VTK file. */
    /* TODO: Migrate to ViennaGrid 3.0
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
            viennashe::write_current_density_to_container(device,
                                                          simulator_obj.potential(),
                                                          quan, viennashe::ELECTRON_TYPE_ID,
                                                          viennashe::models::create_constant_mobility_model(device, 0.1430),
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
            viennashe::write_current_density_to_container(device,
                                                          simulator_obj.potential(),
                                                          quan, viennashe::HOLE_TYPE_ID,
                                                          viennashe::models::create_constant_mobility_model(device, 0.0460),
                                                          current_p);
          }
          my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(current_p), "Hole current density");
        }


        // lattice temperature
        if (quan.get_name() == viennashe::quantity::lattice_temperature())
        {
          typedef typename SimulatorType::SHETimeStepQuantitiesT QuantitiesType;
          typedef typename viennashe::hde::power_density_accessor<DeviceType, QuantitiesType> PowerDensityAccessorType;

          PowerDensityAccessorType pdacc(device,
                                         simulator_obj.quantities(),
                                         simulator_obj.config());

          viennashe::write_macroscopic_quantity_to_container(device, pdacc, pwr_density_container);

          my_vtk_writer.add_scalar_data_on_cells(viennagrid::make_accessor<CellType>(pwr_density_container), "Joule heating power density");
        }


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

    } */

  } //namespace io
} //namespace viennashe

#endif
