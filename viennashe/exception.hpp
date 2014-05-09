#ifndef VIENNASHE_EXCEPTION_HPP
#define VIENNASHE_EXCEPTION_HPP

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

/** @file viennashe/exception.hpp
    @brief Provides the exceptions used in the main viennashe namespace.
*/

namespace viennashe
{

  /** @brief Exception for the case that an invalid cell is encountered */
  class invalid_cell_exception : public std::runtime_error {
  public:
    invalid_cell_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that an invalid doping is encountered */
  class invalid_doping_in_device_exception : public std::runtime_error {
  public:
    invalid_doping_in_device_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that an invalid lattice temperature is encountered */
  class invalid_lattice_temperature_in_device_exception : public std::runtime_error {
  public:
    invalid_lattice_temperature_in_device_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that an invalid lattice temperature is encountered */
  class invalid_trap_energy_exception : public std::runtime_error {
  public:
    invalid_trap_energy_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that an invalid lattice temperature is encountered */
  class invalid_trap_density_exception : public std::runtime_error {
  public:
    invalid_trap_density_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that a segment or a segment that does not belong to a mesh is encountered */
  class invalid_segment_exception : public std::runtime_error {
  public:
    invalid_segment_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that a segment or a segment that does not belong to a mesh is encountered */
  class invalid_boundary_condition_exception : public std::runtime_error {
  public:
    invalid_boundary_condition_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception for the case that an invalid value (depends on the method called) is encountered */
  class invalid_value_exception : public std::runtime_error
  {
  public:
    invalid_value_exception(std::string const & str) : std::runtime_error(str), _value(0) { this->fill_message(); }

    virtual const char* what() const throw()  { return _msg.c_str(); }

    invalid_value_exception(std::string const & str, const double & value)
      : std::runtime_error(str), _value(value)
    { this->fill_message(); }

    virtual ~invalid_value_exception() throw () { }

    private:
      double _value;
      std::string _msg;

      void fill_message()
      {
        std::stringstream ss;
        ss << "Invalid value '" << this->_value << "' encountered. " << std::endl;
        ss << std::runtime_error::what();
        _msg = ss.str();
      }
  };

  /** @brief Exception for the case that a requested feature is not available (due to configuration or due to not having an implementation yet) */
  class unavailable_feature_exception : public std::runtime_error
  {
  public:

    unavailable_feature_exception(std::string const & str) : std::runtime_error(str), _file_name(""), _line("") { this->fill_message(); }

    virtual const char* what() const throw() { return _msg.c_str(); }

    unavailable_feature_exception(std::string const & str, std::string const & file_name, std::string const & line_of_code)
      : std::runtime_error(str), _file_name(file_name), _line(line_of_code)
    { this->fill_message(); }

    virtual ~unavailable_feature_exception() throw () { this->fill_message(); }

    private:
      std::string _file_name;
      std::string _line;
      std::string _msg;

      void fill_message()
      {
        std::stringstream ss;
        ss << "Unavailable feature in '" << _file_name << "' on line '" << _line << "' " << std::endl;
        ss << std::runtime_error::what();
        _msg = ss.str();
      }

  };

  /** @brief Exception in case a (requested) quantity cannot be found */
  class quantity_not_found_exception : public std::runtime_error
  {
  public:
    quantity_not_found_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception thrown in the case that an equation solver failed */
  class solver_failed_exception : public std::runtime_error
  {
  public:
    solver_failed_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception thrown in the case that an equation assembler cannot be found */
  class assembly_not_implemented_exception : public std::runtime_error
  {
  public:
    assembly_not_implemented_exception(std::string const & str) : std::runtime_error(str) {}
  };

  /** @brief Exception thrown in the case that any code does not support the given carrier type */
  class carrier_type_not_supported_exception : public std::runtime_error
  {
  public:
    carrier_type_not_supported_exception(std::string const & str) : std::runtime_error(str) {}
  };


}

#endif



