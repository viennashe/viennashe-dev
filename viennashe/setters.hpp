#ifndef VIENNASHE_SETTERS_HPP
#define VIENNASHE_SETTERS_HPP

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

// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/exception.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/device.hpp"
#include "viennashe/trap_level.hpp"

/** @file  viennashe/setters.hpp
    @brief Contains the definition of convenience functors for accessing device quantities (see class device). Does not contain simulator-specific setters!
*/

namespace viennashe
{

  //
  // Setters
  //


  /** @brief Convenience functor for setting the doping across the device */
  template <typename DeviceType>
  class doping_setter
  {
    public:
      typedef double value_type;

      doping_setter(DeviceType & d, viennashe::carrier_type_id ctype) : device_(d), is_doping_n_(ctype == viennashe::ELECTRON_TYPE_ID) {}

      template <typename T>
      void operator()(T const & t, value_type value) const { is_doping_n_ ? device_.set_doping_n(value, t) : device_.set_doping_p(value, t); }

    private:
      DeviceType & device_;
      bool is_doping_n_;
  };




} // namespace viennashe

#endif

