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

#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/centroid.hpp"


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

    template < typename ElementT >
    double operator()(ElementT const & elem) const
    {
      return this->get_value(elem);
    }

    template < typename ElementT >
    double get_value(ElementT const & elem) const
    {
      if ( values_ != 0)
        return this->values_[elem.id().get()];
      else
        throw viennashe::invalid_value_exception("array_to_accessor: values_ must not be NULL!");
    }

    template < typename ElementT >
    double at(ElementT const & elem) const
    {
      return get_value(elem);
    }


  private:
    double * values_;
  }; // array_to_accessor

  /**
   * @brief The mesh types supported by libviennashe.
   *
   */
  struct meshtype
  {
    enum
    {
      line_1d = 0,
      triangular_2d,
      quadrilateral_2d,
      hexahedral_3d,
      tetrahedral_3d
    };
  };
}


/** @brief Internal C++ to C wrapper for the device. Has typedefs and destructor. */
struct viennashe_device_impl
{
  typedef viennashe::device<viennagrid::line_1d_mesh>          dev1d_type;
  typedef viennashe::device<viennagrid::quadrilateral_2d_mesh> devq2d_type;
  typedef viennashe::device<viennagrid::triangular_2d_mesh>    devt2d_type;
  typedef viennashe::device<viennagrid::hexahedral_3d_mesh>    devh3d_type;
  typedef viennashe::device<viennagrid::tetrahedral_3d_mesh>   devt3d_type;

  viennashe_device_impl() : stype(-1), device_1d(NULL) {  }

  ~viennashe_device_impl()
  {
    if (stype >= 0)
    {
      if (stype == libviennashe::meshtype::line_1d)
        delete device_1d;
      else if (stype == libviennashe::meshtype::quadrilateral_2d)
        delete device_quad_2d;
      else if (stype == libviennashe::meshtype::triangular_2d)
        delete device_tri_2d;
      else if (stype == libviennashe::meshtype::hexahedral_3d)
        delete device_hex_3d;
      else if (stype == libviennashe::meshtype::tetrahedral_3d)
        delete device_tet_3d;

      device_1d = NULL;
    }
    stype = -1;
  }

  bool is_valid() const { return (stype >= 0 && device_1d != NULL); }

  int    stype;

  union
  {
    dev1d_type  * device_1d;
    devq2d_type * device_quad_2d;
    devt2d_type * device_tri_2d;
    devh3d_type * device_hex_3d;
    devt3d_type * device_tet_3d;
  };

}; // viennashe_device_impl


/** @brief Internal C++ to C wrapper for the simulator. Has typedefs and destructor. */
  struct viennashe_simulator_impl
{
  //
  // Mesh typedefs
  //
  typedef viennagrid::line_1d_mesh            mesh1d_type;
  typedef viennagrid::quadrilateral_2d_mesh   meshq2d_type;
  typedef viennagrid::triangular_2d_mesh      mesht2d_type;
  typedef viennagrid::hexahedral_3d_mesh      meshh3d_type;
  typedef viennagrid::tetrahedral_3d_mesh     mesht3d_type;

  //
  // Device typedefs
  //
  typedef viennashe::device<mesh1d_type>  dev1d_type;
  typedef viennashe::device<meshq2d_type> devq2d_type;
  typedef viennashe::device<mesht2d_type> devt2d_type;
  typedef viennashe::device<meshh3d_type> devh3d_type;
  typedef viennashe::device<mesht3d_type> devt3d_type;

  //
  // Simulator typedefs
  //
  typedef viennashe::simulator<dev1d_type>  sim1d_type;
  typedef viennashe::simulator<devq2d_type> simq2d_type;
  typedef viennashe::simulator<devt2d_type> simt2d_type;
  typedef viennashe::simulator<devh3d_type> simh3d_type;
  typedef viennashe::simulator<devt3d_type> simt3d_type;


  viennashe_simulator_impl(int s, viennashe::config * c) : stype(s), conf(c), sim1d(0) {  }

  ~viennashe_simulator_impl()
  {
    if (stype >= 0 && sim1d != NULL)
    {
      if (stype == libviennashe::meshtype::line_1d)
        delete sim1d;
      else if (stype == libviennashe::meshtype::quadrilateral_2d)
        delete simq2d;
      else if (stype == libviennashe::meshtype::triangular_2d)
        delete simt2d;
      else if (stype == libviennashe::meshtype::hexahedral_3d)
        delete simh3d;
      else if (stype == libviennashe::meshtype::tetrahedral_3d)
        delete simt3d;

      conf = NULL; // Not deleted by this wrapper !
      simt3d  = NULL;
    }
    stype = -1;
  }

  bool is_valid() const { return (stype >= 0 && conf != NULL && sim1d != NULL); }

  int     stype;
  viennashe::config * conf; // Not deleted by this wrapper !

  union
  {
    sim1d_type  * sim1d;
    simq2d_type * simq2d;
    simt2d_type * simt2d;
    simh3d_type * simh3d;
    simt3d_type * simt3d;
  };

};

#endif	/* LIBVIENNASHE_VIENNASHE_ALL_HPP */

