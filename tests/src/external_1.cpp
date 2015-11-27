/* =======================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

   Project Head:    Karl Rupp                           rupp@iue.tuwien.ac.at

                    http://viennashe.sourceforge.net/

   License:      MIT (X11), see file LICENSE in the base directory
======================================================================= */

#ifdef _MSC_VER      //Visual Studio complains about potentially dangerous things, which are perfectly legal in our context
  #pragma warning( disable : 4355 )     //use of this in member initializer list
  #pragma warning( disable : 4503 )     //truncated name decoration
#endif

/** \file external_1.cpp Contains checks for the absence of external linkage (otherwise, library is not truly 'header-only')
 *  \test A check for the absence of external linkage (otherwise, library is not truly 'header-only')
 */



//
// *** System
//
#include <iostream>

//#define VIENNASHE_HAVE_PARALLEL_SOLVER    //automatically set by the build system
//#define VIENNASHE_HAVE_GPU_SOLVER         //automatically set by the build system

#include "viennagrid/viennagrid.h"


#include "viennashe/core.hpp"


//defined in external_2.cpp
void other_func();

//
// -------------------------------------------------------------
//
int main()
{
  std::cout << "*****************" << std::endl;
  std::cout << "* Test started! *" << std::endl;
  std::cout << "*****************" << std::endl;

  //this is the external linkage check:
  other_func();

  std::cout << "*******************************" << std::endl;
  std::cout << "* Test finished successfully! *" << std::endl;
  std::cout << "*******************************" << std::endl;

  return EXIT_SUCCESS;
}
