#ifndef VIENNASHE_IO_INITIALGUESS_WRITER_HPP
#define VIENNASHE_IO_INITIALGUESS_WRITER_HPP

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


/** @file viennashe/io/initial_guess_writer.hpp
    @brief A writer function for writing the initial guess to an output file
*/

#include <vector>
#include "viennagrid/viennagrid.h"

namespace viennashe
{
  namespace io
  {

    /**
     * @brief Writes the inital guess to a VTK file
     * @param mesh    The mesh on which the guess is defined
     * @param key     The key to the inital guess used in the VTK file
     * @param qas     The quantity accessor
     */
    template<typename MeshT, typename QuanAssembler>
    void write_initial_guess(MeshT const& mesh, std::string const& key, QuanAssembler qas)
    {

      viennagrid_element_id *vertices_begin, *vertices_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, 0, &vertices_begin, &vertices_end));

      std::vector<double> vtk_data(vertices_end - vertices_begin);
      for(viennagrid_element_id *vit = vertices_begin; vit != vertices_end; vit++)
      {
        vtk_data[viennagrid_index_from_element_id(*vit)] = qas.current_guess_at(*vit);
      }

      //viennagrid::io::vtk_writer<MeshT> my_vtk_writer;
      //my_vtk_writer.add_scalar_data_on_vertices(viennagrid::make_accessor<VertexType>(vtk_data), key);
      //my_vtk_writer(mesh, key);
    }

  } // io
} // viennashe


#endif

