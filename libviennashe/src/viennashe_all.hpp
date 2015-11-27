#ifndef LIBVIENNASHE_VIENNASHE_ALL_HPP
#define	LIBVIENNASHE_VIENNASHE_ALL_HPP
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

/** @file libviennashe/src/viennashe_all.hpp
    @brief Contains all viennashe includes and most of the C++/C-wrappers. Contains macros
*/

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include "libviennashe/include/sys.h"

// ++++++++++++++++++++++++++++++++++++
//
// +++++++++++   C++ ONLY   +++++++++++
//
// ++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cstdlib>
#include <vector>

#include "viennashe/forwards.h"
#include "viennashe/device.hpp"
#include "viennashe/simulator.hpp"
#include "viennashe/io/she_vtk_writer.hpp"
#include "viennashe/io/vector.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/config.hpp"
#include "viennashe/she/postproc/all.hpp"
#include "viennashe/io/gnuplot_writer.hpp"

#include "viennashe/postproc/current_density.hpp"
#include "viennashe/postproc/electric_field.hpp"
#include "viennashe/postproc/electric_flux_density.hpp"

#include "viennashe/models/all.hpp"

#include "viennagrid/viennagrid.h"


#ifndef __FUNCTION_NAME__
    #ifdef WIN32   // WINDOWS
        #define __FUNCTION_NAME__   __FUNCTION__
    #else          // UNIX
        #define __FUNCTION_NAME__   __func__
    #endif
#endif


/* ============================================================================================================ */
#define CHECK_ARGUMENT_FOR_NULL(arg, pos, name) if((arg) == NULL) \
                                                { \
                                                  viennashe::log::error() << "ERROR: " << __FUNCTION_NAME__ \
                                                  << "(): " << (name) << " must not be NULL!" << std::endl; \
                                                  return (pos); \
                                                }
/* ============================================================================================================ */
#define CHECK_ARGUMENT(mycheck, pos, msg) if(mycheck) \
                                           { \
                                             viennashe::log::error() << "ERROR: " << __FUNCTION_NAME__ \
                                             << "(): " << (msg) << "!" << std::endl; \
                                             return (pos); \
                                           }
/* ============================================================================================================ */

struct viennashe_device_impl
{
  viennashe::device<viennagrid_mesh> device_;
};

struct viennashe_simulator_impl
{
  viennashe_simulator_impl(viennashe::device<viennagrid_mesh> dev, viennashe::config conf) : sim_(dev, conf) {}

  viennashe::simulator<viennashe::device<viennagrid_mesh> > sim_;
};

/**
 * @brief The internal C++ namespace of the library
 */
namespace libviennashe
{
  /** @brief Maps a C-array to an element based accessor. Uses element.id() as an index to the C-array */
  class array_to_accessor
  {
  public:

    array_to_accessor(double * vals) : values_(vals) { }

    double operator()(viennagrid_element_id elem) const { return get_value(elem); }
    double         at(viennagrid_element_id elem) const { return get_value(elem); }

    double get_value(viennagrid_element_id elem) const
    {
      if (values_)
        return this->values_[viennagrid_index_from_element_id(elem)];
      else
        throw viennashe::invalid_value_exception("array_to_accessor: values_ must not be NULL!");
    }

  private:
    double * values_;
  }; // array_to_accessor

}


#endif	/* LIBVIENNASHE_VIENNASHE_ALL_HPP */

