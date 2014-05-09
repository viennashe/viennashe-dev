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
#include "viennagrid/accessor.hpp"
#include "viennagrid/io/vtk_writer.hpp"

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
      typedef typename viennagrid::result_of::vertex<MeshT>::type                   VertexType;
      typedef typename viennagrid::result_of::const_vertex_range<MeshT>::type       VertexContainer;
      typedef typename viennagrid::result_of::iterator<VertexContainer>::type       VertexIterator;

      VertexContainer vertices(mesh);
      std::vector<double> vtk_data(vertices.size());
      for(VertexIterator vit = vertices.begin(); vit != vertices.end(); vit++)
      {
        vtk_data[vit->id().get()] = qas.current_guess_at(*vit);
      }

      viennagrid::io::vtk_writer<MeshT> my_vtk_writer;
      my_vtk_writer.add_scalar_data_on_vertices(viennagrid::make_accessor<VertexType>(vtk_data), key);
      my_vtk_writer(mesh, key);
    }

  } // io
} // viennashe


#endif

