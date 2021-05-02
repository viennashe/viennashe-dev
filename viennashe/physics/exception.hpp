#ifndef VIENNASHE_PHYSICS_EXCEPTION_HPP
#define VIENNASHE_PHYSICS_EXCEPTION_HPP

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

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

/** @file viennashe/physics/exception.hpp
    @brief All the exceptions used within the viennashe::physics namespace.
*/

namespace viennashe
{
  namespace physics
  {

    /** @brief Exception which is thrown if there is no lattice temperature for a vertex */
    class no_lattice_temperature_available_exception : public std::runtime_error
    {
    public:

      virtual const char* what() const throw() { return msg_.c_str(); }

      no_lattice_temperature_available_exception(std::string const & str,
                                                 const double & temp,
                                                 const long & vertex_id)
        : std::runtime_error(str), temp_(temp), vt_id_(vertex_id)  { this->fill_message(); }

      virtual ~no_lattice_temperature_available_exception() throw() { }

    private:
      double temp_;
      long vt_id_;

      std::string msg_;

      void fill_message()
      {
        std::stringstream ss;
        ss << "Invalid temperature '" << temp_ << "' on vertex '" << vt_id_ << "'" << std::endl;
        msg_ = ss.str();
      }

    };


    /** @brief Exception which is thrown if there is no data for a given energy */
    class data_not_available_exception : public std::runtime_error
    {
    public:

      virtual const char* what() const throw() { return _msg.c_str(); }

      data_not_available_exception(std::string const & str,
                                   const std::size_t & bin_no)
        : std::runtime_error(str), _bin_no(bin_no)
      { }

      virtual ~data_not_available_exception() throw() { }

    private:
      std::size_t _bin_no;
      std::string _msg;

      void fill_message()
      {
        std::stringstream ss;
        ss << std::runtime_error::what() << " ## No data found for data bin '" << this->_bin_no << "'" << std::endl;
        _msg = ss.str();
      }
    };


    /** @brief Exception which is thrown if there is no data for a given energy */
    class dispersion_not_invertible_exception : public std::runtime_error
    {
    public:
      dispersion_not_invertible_exception(std::string const & str) : std::runtime_error(str) { }
    };

  } // namespace io
} // namespace viennashe


#endif

