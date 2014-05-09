#ifndef VIENNASHE_IO_GNUPLOT_WRITER_HPP
#define VIENNASHE_IO_GNUPLOT_WRITER_HPP

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
#include "viennashe/physics/constants.hpp"
#include "viennashe/io/exception.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/log/log.hpp"


/** @file viennashe/io/gnuplot_writer_edf.hpp
    @brief Writes the energy distribution function to a file which can be processed by Gnuplot.
 */

namespace viennashe
{
  namespace io
  {

    /** @brief Writes quantities to a file which can be processed by Gnuplot. Works for 1d, 2d and 3d only. */
    struct gnuplot_writer
    {

      /** @brief Helper function for writing a quantity at particular cells to a plain text file. Output can be processed by e.g. Gnuplot */
      template <typename DeviceType,
                typename CellFilterType,
                typename QuantitiyAccessorType>
      void operator()(DeviceType const & device,
                      CellFilterType const & cell_filter,
                      QuantitiyAccessorType const & quan,
                      const std::string filename) const
      {
        typedef typename DeviceType::mesh_type  MeshType;

        typedef typename viennagrid::result_of::point<MeshType>::type     PointType;

        typedef typename viennagrid::result_of::const_cell_range<MeshType>::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type       CellIterator;

        std::ofstream writer(filename.c_str());

        if ( !writer )
        {
          throw viennashe::io::cannot_open_file_exception(filename);
        }

        writer << "## ViennaSHE - gnuplot output " << std::endl;

        CellContainer cells(device.mesh());
        for ( CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit )
        {
          // filter vertices
          if ( !cell_filter(*cit) ) continue;

          //write values at point
          PointType p = viennagrid::centroid(*cit);
          for (std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i) writer << p[i] << " ";
          this->write_value(writer, quan(*cit));
          writer << std::endl;
        } // for vertices

        log::info() << "* write_vertex_quantity_for_gnuplot(): Writing data to '"
                  << filename
                  << "' (can be viewed with e.g. gnuplot)" << std::endl;
      } // operator()

      /** @brief Helper function for writing a quantity at all cells to a plain text file. Output can be processed by e.g. Gnuplot */
      template <typename DeviceType,
                typename QuantitiyAccessorType>
      void operator()(DeviceType const & device,
                      QuantitiyAccessorType const & quan,
                      const std::string filename) const
      {
        this->operator ()(device, viennashe::util::any_filter(), quan, filename);
      }


    private:

      template < typename WriterT, typename ValueT >
      void write_value(WriterT &, ValueT const &) const
      {
        typedef typename error_indicator<ValueT>::ERROR_NO_IMPLEMENTATION_FOR_THE_GIVEN_VALUET   error_type;
        error_type * obj; obj = 0; // fixes compiler warning; does not provide any functionality
      }

      template < typename WriterT >
      void write_value(WriterT & writer, double const & val) const
      {
        writer << val << " " ;
      }

      template < typename WriterT >
      void write_value(WriterT & writer, std::vector<double> const & val) const
      {
        for (std::size_t i = 0; i < val.size(); ++i)
          writer << val[i] << " ";
      }


    }; // gnuplot_writer

    /** @brief Writes a quantity (on vertices) per point to a text file suitable for gnuplot */
    template < typename DeviceType, typename AccessorType >
    void write_cell_quantity_for_gnuplot(AccessorType const & quan, DeviceType const & device, std::string filename)
    {
      gnuplot_writer my_writer;
      my_writer(device, viennashe::util::any_filter(), quan, filename);
    } // write_vertex_quantity_for_gnuplot())

  } //namespace io

} //namespace viennashe

#endif
