#ifndef VIENNASHE_IO_GNUPLOT_WRITER_EDF_HPP
#define VIENNASHE_IO_GNUPLOT_WRITER_EDF_HPP

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
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/centroid.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/config.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/io/exception.hpp"

/** @file viennashe/io/gnuplot_writer_edf.hpp
    @brief Writes the energy distribution function to a file which can be processed by Gnuplot.
 */

namespace viennashe
{
  namespace io
  {

    /** @brief Writes the energy distribution function to a file which can be processed by Gnuplot. Works for 1d, 2d and 3d only. */
    struct gnuplot_edf_writer
    {

      /** @brief Triggers the write process
       *
       * @param device         The device (includes a ViennaGrid mesh) on which simulation is carried out
       * @param edf            The distribution function
       * @param coordinate_x   x-coordinate of the point for which the data should be written
       * @param filename       Name of the file to be written to
       */
      template < typename DeviceType, typename EDFWrapperT >
      void operator()(DeviceType const & device,
                      EDFWrapperT const & edf,
                      const double coordinate_x,
                      const std::string filename) const
      {
        typedef typename DeviceType::mesh_type            MeshType;
        typedef typename viennagrid::result_of::const_cell_range<MeshType> ::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type        CellIterator;
        typedef typename EDFWrapperT::she_quantity_type         she_quantity_type;

        she_quantity_type        const & quan       = edf.quan();
        //dispersion_relation_type const & dispersion = edf.dispersion_relation();

        std::map<std::pair<double, double>, double > values; //key is x-coordinate, value is velocity
        std::ofstream writer(filename.c_str());

        if ( !writer )
          throw cannot_open_file_exception(filename);

        //iterate over edges:
        CellContainer cells(device.mesh());
        for ( CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit )
        {
          if ( viennagrid::centroid(*cit)[0] != coordinate_x )
            continue;

          for ( std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H )
          {
            const long index = quan.get_unknown_index(*cit, index_H);
            const double eps = quan.get_kinetic_energy(*cit, index_H);
            if ( index > -1 )
              values[std::make_pair(viennagrid::centroid(*cit)[1], eps)] = edf(*cit, eps, index_H);
          }
        }


        //now stream values to file:
        double last_x = -1.0;
        for ( std::map<std::pair<double, double>, double >::iterator it = values.begin();
              it != values.end();
              ++it )
        {
          double x_coord = it->first.first;
          if ( x_coord < last_x || x_coord > last_x )
          {
            writer << std::endl;
            last_x = x_coord;
          }
          writer << it->first.first << " " << it->first.second << " " << it->second << std::endl;
        }

        writer.close();
      }

      /** @brief Helper function for writing the energy distribution function at a particular points to a plain text file. Output can be processed by e.g. Gnuplot
       * @param device The device
       * @param cell_filter A cell filter (has operator(CellType const & cell)), which returns true if the writer should be doing something for the given cell
       * @param edf The energy distribution function (cf. df_wrappers.hpp)
       * @param filename The path/filename where the data gets written to
       */
      template <typename DeviceType,
                typename CellFilterType,
                typename EDFWrapperT>
      void operator()(DeviceType const & device,
                      CellFilterType const & cell_filter,
                      EDFWrapperT const & edf,
                      const std::string filename)
      {
        typedef typename DeviceType::mesh_type   MeshType;

        typedef typename viennagrid::result_of::point<MeshType>::type   PointType;
        typedef typename viennagrid::result_of::const_cell_range<MeshType> ::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type        CellIterator;

        typedef typename EDFWrapperT::dispersion_relation_type  dispersion_relation_type;
        typedef typename EDFWrapperT::she_quantity_type         she_quantity_type;

        she_quantity_type        const & quan       = edf.quan();
        dispersion_relation_type const & dispersion = edf.dispersion_relation();

        std::ofstream writer(filename.c_str());

        if ( !writer )
        {
          throw viennashe::io::cannot_open_file_exception(filename);
        }

        CellContainer cells(device.mesh());
        for ( CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit )
        {
          // filter vertices
          if ( !cell_filter(*cit) ) continue;

          //write values at point
          for ( std::size_t index_H = 1;
                index_H < quan.get_value_H_size() - 1;
                ++index_H )
          {

            const long index = quan.get_unknown_index(*cit, index_H);
            const double energy = quan.get_kinetic_energy(*cit, index_H);
            if ( index > -1 && energy > 0 )
            {
              for (std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i) writer << viennagrid::centroid(*cit)[i] << " ";
              writer << (energy / viennashe::physics::constants::q) << " "
                     << edf(*cit, energy, index_H) << " "
                     << dispersion.density_of_states(energy)
                     << std::endl;
            }
          } // for index_H
          writer << std::endl << std::endl; // for gnuplot block index
        } // for vertices
      } // operator()

    };

  } //namespace io
} //namespace viennashe

#endif
