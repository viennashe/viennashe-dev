#ifndef VIENNASHE_IO_EXCEPTION_HPP_
#define VIENNASHE_IO_EXCEPTION_HPP_

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


#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

/** @file viennashe/io/exception.hpp
    @brief All the exceptions used within the viennashe::io namespace.
*/

namespace viennashe
{
  /** @brief Namespace for all (file-)IO specific code  */
  namespace io
  {

    //
    // File IO-specific
    //
    /** @brief Exception which is thrown if a file cannot be opened */
    class cannot_open_file_exception : public std::runtime_error {
    public:
      cannot_open_file_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown if a End-of-File is encountered even though further data is expected */
    class premature_end_of_file_exception : public std::runtime_error {
    public:
      premature_end_of_file_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown if a dimension different from the current dimension of viennagrid is found in the file.  */
    class invalid_dimension_exception : public std::runtime_error {
    public:
      invalid_dimension_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown whenever a requested IO feature is currently not supported.  */
    class io_operation_unsupported_exception : public std::runtime_error {
    public:
      io_operation_unsupported_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown whenever a requested attribute was not found in a file.  */
    class file_attribute_not_found_exception : public std::runtime_error {
    public:
      file_attribute_not_found_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception which is thrown whenever a requested segment was not found in a file.  */
    class file_segment_not_found_exception : public std::runtime_error {
    public:
      file_segment_not_found_exception(std::string const & str) : std::runtime_error(str) {}
    };

  } //namespace io
}//namespace viennashe



#endif
