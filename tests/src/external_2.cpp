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

#ifdef _MSC_VER      //Visual Studio complains about potentially dangerous things, which are perfectly legal in our context
  #pragma warning( disable : 4355 )     //use of this in member initializer list
  #pragma warning( disable : 4503 )     //truncated name decoration
#endif

/** \file external_2.cpp Yet another check of external linkage (we require "header-only" capability)
 *  \test A second check for the absence of external linkage (otherwise, library is not truly 'header-only')
 */

#include <iostream>

//#define VIENNASHE_HAVE_PARALLEL_SOLVER    //automatically set by the build system
//#define VIENNASHE_HAVE_GPU_SOLVER         //automatically set by the build system

#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/config/default_configs.hpp"

#include "viennashe/core.hpp"

void other_func();

void other_func()
{
  viennagrid::tetrahedral_3d_mesh  mesh;

  std::cout << "--- Tetrahedral mesh, 3d ---" << std::endl;
  std::cout << "Size<0>: " << viennagrid::vertices(mesh).size() << std::endl;
  std::cout << "Size<1>: " << viennagrid::vertices(mesh).size() << std::endl;
  std::cout << "Size<2>: " << viennagrid::vertices(mesh).size() << std::endl;
  std::cout << "Size<3>: " << viennagrid::vertices(mesh).size() << std::endl;
}
