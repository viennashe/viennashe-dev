#ifndef VIENNASHE_IO_VECTOR_HPP_
#define VIENNASHE_IO_VECTOR_HPP_

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

#include "viennashe/io/exception.hpp"

/** @file viennashe/io/vector.hpp
    @brief Simple routines for reading a vector from file, or writing a vector to file.
*/

namespace viennashe
{
  namespace io
  {

    /** @brief Reads a vector from a file
     *
     * @param vec      The vector
     * @param filename The filename
     */
    template <typename VectorType>
    void read_vector_from_file(VectorType & vec,
                               const std::string & filename)
    {
      std::ifstream file(filename.c_str());

      if (!file)
        throw cannot_open_file_exception(filename);

      std::size_t size;
      file >> size;

      //if(size == 0) throw std::range_error(std::string("There are zero vector-components in File: '") + filename + std::string("'"))

      vec.resize(size, false);

      for (std::size_t i = 0; i < size; ++i)
      {
        double element;
        if (file.eof())
          throw premature_end_of_file_exception(filename);
        file >> element;
        vec[i] = element;
      }
    }

    /** @brief Writes a vector to a file
     *
     * @param vec      The vector
     * @param filename The filename
     */
    template <typename VectorType>
    void write_vector_to_file(VectorType & vec,
                              const std::string & filename)
    {
      std::ofstream writer(filename.c_str());

      if (!writer)
        throw cannot_open_file_exception(filename);

      writer << vec.size() << std::endl;
      for (std::size_t i = 0; i < vec.size(); ++i)
      {
        writer << vec[i] << std::endl;
      }
    }


  } //namespace io
}//namespace viennashe



#endif
