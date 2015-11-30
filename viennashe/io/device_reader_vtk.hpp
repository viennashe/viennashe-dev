#ifndef VIENNASHE_IO_DEVICE_READER_VTK_HPP
#define VIENNASHE_IO_DEVICE_READER_VTK_HPP

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

// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/log/log.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/io/exception.hpp"

#include "viennashe/materials/all.hpp"


namespace viennashe
{
  namespace io
  {

    namespace detail
    {

      struct mesh_generator_vtk
      {
        mesh_generator_vtk(std::string filename) : filename_(filename)
        {
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_create(&mesh_reader_));
        }

        ~mesh_generator_vtk() { VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_release(mesh_reader_)); }

        void operator()(viennagrid_mesh mesh)
        {
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_mesh_set(mesh_reader_, mesh));
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_read(mesh_reader_, filename_.c_str()));
        }

        viennagrid_mesh_io & reader() { return this->mesh_reader_; }

        private:
          viennagrid_mesh_io mesh_reader_;
          std::string filename_;
      };

    }

    /**
     * @brief Reads and initalises a device from a VTK file. If material_key is empty the material IDs will not be initalised!
     *
     * @param device The device to initalise
     * @param filename The name/path to the VTK file
     * @param doping_n_key The VTK key to the donor doping
     * @param doping_p_key The VTK key to the acceptor doping
     * @param material_key The VTK key to the material id
     */
    template < typename DeviceType >
    bool read_device_vtk(DeviceType & device,
                         std::string filename,
                         std::string doping_n_key,
                         std::string doping_p_key,
                         std::string material_key)
    {
      if (doping_n_key.empty() || doping_p_key.empty())
      {
        viennashe::log::error() << "* read_device(): One of the doping keys is empty! Unable to proceed." << std::endl;
        return false;
      }

      detail::mesh_generator_vtk my_vtk_reader(filename);

      device.load_device(my_vtk_reader); // basic init of mesh and device

      /* TODO: Migrate code

      viennagrid_mesh_io & mesh_reader = my_vtk_reader.reader();

      const std::size_t num_segments = device.segmentation().size();
      log::info() << "* read_device_vtk(): There are " << num_segments << " segments." << std::endl;

      // DEBUG INFORMATION:
      for (std::size_t j = 0; j < num_segments; ++j)
      {
        log::debug() << "* read_device():" << std::endl;

        if(mesh_reader.scalar_cell_data_names(j).size() > 0)
          log::debug() << "  Segment " << j << " has the following scalar CELL quantities: " << std::endl;
        for(std::size_t i = 0; i < mesh_reader.scalar_cell_data_names(j).size(); ++i)
          log::debug() << "\t" << mesh_reader.scalar_cell_data_names(j)[i] << std::endl;

        if(mesh_reader.scalar_vertex_data_names(j).size() > 0)
          log::debug() << "  Segment " << j << " has the following scalar VERTEX quantities: " << std::endl;
        for(std::size_t i = 0; i < mesh_reader.scalar_vertex_data_names(j).size(); ++i)
          log::debug() << "\t" << mesh_reader.scalar_vertex_data_names(j)[i] << std::endl;

        log::debug() << std::endl;
      }

      //
      // Set the material per cell on each segment
      if (!material_key.empty())
      {
        viennashe::log::info() << "* read_device(): Setting materials ..." << std::endl;
        log::debug()<< "* read_device(): Using '" << material_key << "' as material ID per cell." << std::endl;

        for (std::size_t j = 0; j < num_segments; ++j)
        {
          bool has_material_information = false;
          std::vector<std::string> cell_data_names = mesh_reader.scalar_cell_data_names(j);
          for (std::size_t n=0; n<cell_data_names.size(); ++n)
          {
            if (cell_data_names[n] == material_key)
              has_material_information = true;
          }

          if (!has_material_information)
          {
            viennashe::log::error() << "* read_device(): WARNING: Unable to find material for segment " << j << "! Skipping..." << std::endl;
            continue;
          }

          scalar_cell_data const & mat_data = mesh_reader.cell_scalar_field(material_key, j);
          CellContainer cells(device.segment(j));
          for ( CellIterator cit = cells.begin(); cit != cells.end(); ++cit )
          {
            const long material_id = mat_data.at(*cit);
            device.set_material(material_id, *cit);
          } // for cells
        }
      }
      else
      {
        viennashe::log::warning() << "* read_device(): Empty material key found! Assuming the material is being set by the user ... " << std::endl;
      }


      viennashe::log::info() << "* read_device(): Setting doping ..." << std::endl;

      for (std::size_t j = 0; j < num_segments; ++j)
      {
        bool doping_on_cells = false;

        std::vector<std::string> cell_data_names = mesh_reader.scalar_cell_data_names(j);
        for (std::size_t n=0; n<cell_data_names.size(); ++n)
        {
          if (cell_data_names[n] == doping_n_key)
            doping_on_cells = true;
        }

        if (doping_on_cells && mesh_reader.cell_scalar_field(doping_n_key, j).is_valid() ) // look for doping on cells
        {
          log::debug()<< "* read_device(): Using '" << doping_n_key << "' as donor doping per cell on segment " << j << "." << std::endl;
          log::debug()<< "* read_device(): Using '" << doping_p_key << "' as acceptor doping per cell on segment " << j << "." << std::endl;

          scalar_cell_data const & data_n = mesh_reader.cell_scalar_field(doping_n_key, j);
          scalar_cell_data const & data_p = mesh_reader.cell_scalar_field(doping_p_key, j);

          if (!data_n.is_valid() || !data_p.is_valid())
          {
            log::error() << "* read_device(): Invalid doping key on cells !" << std::endl;
            return false;
          }

          CellContainer cells(device.segment(j));
          for ( CellIterator cit = cells.begin(); cit != cells.end(); ++cit )
          {
            device.set_doping_n(data_n.at(*cit), *cit);
            device.set_doping_p(data_p.at(*cit), *cit);
          } // for cells
        }
        else if (!doping_on_cells && mesh_reader.vertex_scalar_field(doping_n_key, j).is_valid() ) // doping on vertices
        {
          log::debug()<< "* read_device(): Using '" << doping_n_key << "' as donor doping per vertex." << std::endl;
          log::debug()<< "* read_device(): Using '" << doping_p_key << "' as acceptor doping per vertex." << std::endl;

          scalar_vertex_data const & data_n = mesh_reader.vertex_scalar_field(doping_n_key, j);
          scalar_vertex_data const & data_p = mesh_reader.vertex_scalar_field(doping_p_key, j);

          if (!data_n.is_valid() || !data_p.is_valid())
          {
            log::error() << "* read_device(): Invalid doping key on vertices !" << std::endl;
            return false;
          }

          CellContainer cells(device.segment(j));
          for ( CellIterator cit = cells.begin(); cit != cells.end(); ++cit )
          {
            double doping_n = 1.0;
            double doping_p = 1.0;
            VertexOnCellContainer vertices(*cit);
            double N_n = 0;
            double N_p = 0;
            for ( VertexOnCellIterator vocit = vertices.begin(); vocit != vertices.end(); ++vocit )
            {
              if (data_n.at(*vocit))
              {
                doping_n *= data_n.at(*vocit);
                ++N_n;
              }
              if (data_p.at(*vocit))
              {
                doping_p *= data_p.at(*vocit);
                ++N_p;
              }
            }
            if (N_n > 0)
              device.set_doping_n(std::pow(doping_n, 1.0/N_n), *cit);
            if (N_p > 0)
              device.set_doping_p(std::pow(doping_p, 1.0/N_p), *cit);
          } // for cells
        }
        else
        {
          viennashe::log::error() << "* read_device(): None or inconsistent (cells,vertices) donor doping found! " << std::endl;
          return false;
        }
      } */

      return true;
    } // read_device

    template <typename DeviceT, typename SimulatorT>
    bool read_initial_guess_vtk(DeviceT const & device,
                                SimulatorT & simulator,
                                std::string filename,
                                std::string quantity_name,
                                std::string vtk_quantity_key)
    {
      // TODO: migrate code
      /*
      detail::mesh_generator_vtk<MeshType> my_vtk_reader(filename);

      typename detail::mesh_generator_vtk<MeshType>::vtk_reader_type & mesh_reader = my_vtk_reader.reader();

      DeviceT dummy_device;
      dummy_device.load_device(my_vtk_reader); // basic init of mesh and device

      const std::size_t num_segments = dummy_device.segmentation().size();
      log::info() << "* read_initial_guess_vtk(): There are " << num_segments << " segments." << std::endl;

      // DEBUG INFORMATION:
      for (std::size_t j = 0; j < num_segments; ++j)
      {
        log::debug() << "* read_initial_guess_vtk():" << std::endl;

        if(mesh_reader.scalar_cell_data_names(j).size() > 0)
          log::debug() << "  Segment " << j << " has the following scalar CELL quantities: " << std::endl;
        for(std::size_t i = 0; i < mesh_reader.scalar_cell_data_names(j).size(); ++i)
          log::debug() << "\t" << mesh_reader.scalar_cell_data_names(j)[i] << std::endl;

        if(mesh_reader.scalar_vertex_data_names(j).size() > 0)
          log::debug() << "  Segment " << j << " has the following scalar VERTEX quantities: " << std::endl;
        for(std::size_t i = 0; i < mesh_reader.scalar_vertex_data_names(j).size(); ++i)
          log::debug() << "\t" << mesh_reader.scalar_vertex_data_names(j)[i] << std::endl;

        log::debug() << std::endl;
      }

      viennashe::log::info() << "* read_initial_guess_vtk(): Setting quantity as initial guess ..." << std::endl;

      std::deque<double> cell_data(viennagrid::cells(device.mesh()).size());
      for (std::size_t j = 0; j < num_segments; ++j)
      {
        bool quantity_on_cells = false;

        std::vector<std::string> cell_data_names = mesh_reader.scalar_cell_data_names(j);
        for (std::size_t n=0; n<cell_data_names.size(); ++n)
        {
          if (cell_data_names[n] == vtk_quantity_key)
            quantity_on_cells = true;
        }

        if (quantity_on_cells && mesh_reader.cell_scalar_field(vtk_quantity_key, j).is_valid() ) // look for quantity on cells
        {
          log::debug()<< "* read_initial_guess_vtk(): Using '" << vtk_quantity_key << "' as VTK quantity key for " << quantity_name << " per cell on segment " << j << "." << std::endl;

          scalar_cell_data const & data = mesh_reader.cell_scalar_field(vtk_quantity_key, j);

          if (!data.is_valid())
          {
            log::error() << "* read_initial_guess_vtk(): Invalid quantity key on cells !" << std::endl;
            return false;
          }

          CellContainer cells(dummy_device.segment(j));
          for ( CellIterator cit = cells.begin(); cit != cells.end(); ++cit )
            cell_data.at(cit->id().get()) = data.at(*cit);
        }
        else if (!quantity_on_cells && mesh_reader.vertex_scalar_field(vtk_quantity_key, j).is_valid() ) // doping on vertices
        {
          log::debug()<< "* read_initial_guess_vtk(): Using '" << vtk_quantity_key << "' as VTK quantity key for " << quantity_name << " per vertex on segment " << j << "." << std::endl;

          scalar_vertex_data const & data = mesh_reader.vertex_scalar_field(vtk_quantity_key, j);

          if (!data.is_valid())
          {
            log::error() << "* read_device(): Invalid quantity key on vertices !" << std::endl;
            return false;
          }

          CellContainer cells(dummy_device.segment(j));
          for ( CellIterator cit = cells.begin(); cit != cells.end(); ++cit )
          {
            double value = 0.0;
            VertexOnCellContainer vertices(*cit);
            const double N = static_cast<double>( vertices.size() );
            for ( VertexOnCellIterator vocit = vertices.begin(); vocit != vertices.end(); ++vocit )
              value += data.at(*vocit);
            cell_data.at(cit->id().get()) = value / N;
          } // for cells
        }
        else
        {
          viennashe::log::error() << "* read_initial_guess_vtk(): None or inconsistent (cells,vertices) quantity data for " << quantity_name << " found! " << std::endl;
          return false;
        }
      }

      typename viennagrid::result_of::accessor<std::deque<double>, CellType>::type cell_data_wrapper(cell_data);
      simulator.set_initial_guess(quantity_name, cell_data_wrapper);
      */

      return true;
    } // read_device

    /**
     * @brief Reads and initalises a device from a VTK file. Assumes that the doping and material can be found using "doping_n", "doping_p", "material".
     *
     * @param device The device to initalise
     * @param filename The name/path to the VTK file
     */
    template < typename DeviceType >
    void read_device_vtk(DeviceType & device, const std::string filename)
    {
      read_device(device, filename, "doping_n", "doping_p", "material");
    } // read_device

  } // io
} // viennashe

#endif /* VIENNASHE_IO_DEVICE_READER_VTK_HPP */

