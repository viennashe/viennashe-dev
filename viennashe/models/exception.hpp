#ifndef VIENNASHE_MODELS_EXCEPTION_HPP
#define VIENNASHE_MODELS_EXCEPTION_HPP
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


#include <stdexcept>
#include <sstream>

/** @file viennashe/models/exception.hpp
    @brief Contains exceptions for viennashe::models
 */


namespace viennashe
{
  namespace models
  {
    /** @brief Thrown whenever an invalid parameter to a model is given */
    class invalid_parameter_exception : public std::runtime_error
    {
      public:
        invalid_parameter_exception(std::string const & str)
            : std::runtime_error(str), msg_(""), value_(0) { this->fill_message(); }

        virtual const char* what() const throw() { return this->msgcache_.c_str();  }

        invalid_parameter_exception(std::string msg, double value)
            : std::runtime_error(msg), msg_(msg), value_(value)  { this->fill_message(); }

        virtual ~invalid_parameter_exception() throw() {}

    private:

      void fill_message()
      {
        std::stringstream ss;
        ss << "* Models: Invalid parameter! '" << msg_ << "' Invalid value being: '" << value_ << "'";
        msgcache_ = ss.str();
      }

      std::string msg_;
      double value_;
      std::string msgcache_;
    };

    /** @brief Thrown whenever a model finds an index to be out of bounds */
    class index_out_of_bounds_exception : public std::runtime_error
    {
      public:
        index_out_of_bounds_exception(std::string const & str)
            : std::runtime_error(str), msg_(""), index_(0), bound_(0), msgcache_("") { }

        virtual const char* what() const throw() { return this->msgcache_.c_str(); }

        index_out_of_bounds_exception(std::string msg, long index, long bound)
            : std::runtime_error(msg), msg_(msg), index_(index), bound_(bound)
         { this->fill_message(); }

        virtual ~index_out_of_bounds_exception() throw() {}

    private:

      void fill_message()
      {
        std::stringstream ss;
        ss << "* Models: index out of bounds! '" << msg_ << "' index being: '" << index_ << "'" << "' bound being: '" << bound_ << "'";
        this->msgcache_ = ss.str();
      }

      std::string msg_;
      long index_;
      long bound_;
      std::string msgcache_;
    };

    /** @brief Exception for the case that the evaluation of a model fails */
    class model_evaluation_exception : public std::runtime_error {
    public:
      model_evaluation_exception(std::string const & str) : std::runtime_error(str) {}
    };

  } // namespace models

} // namespace viennashe

#endif /* VIENNASHE_MODELS_EXCEPTION_HPP */

