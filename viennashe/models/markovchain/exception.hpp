#ifndef VIENNASHE_MODELS_MARKOVCHAIN_EXCEPTION_HPP
#define VIENNASHE_MODELS_MARKOVCHAIN_EXCEPTION_HPP
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
#include "viennashe/models/exception.hpp"


/** @file viennashe/models/markovchain/exception.hpp
    @brief Contains exceptions for markov-models in viennashe::models
 */

namespace viennashe
{
  namespace models
  {

    /** @brief Thrown whenever a markov chain model finds an invalid or non-existing state */
    class invalid_state_exception : public std::runtime_error
    {
      public:
        invalid_state_exception(std::string const & str)
            : std::runtime_error(str), msg_(""), index_(0), msgcache_("") { }

        virtual const char* what() const throw() { return this->msgcache_.c_str(); }

        invalid_state_exception(std::string msg, std::size_t index)
            : std::runtime_error(msg), msg_(msg), index_(index)
         { this->fill_message(); }

        virtual ~invalid_state_exception() throw() {}

    private:

      void fill_message()
      {
        std::stringstream ss;
        ss << "* Models: invalid state '" << index_ << "' does not exist or is malconfigured, i.e. state.occupancy() > 1.0 or state.occupancy() < 0 !";
        ss << "  Message: '" << msg_ << "'";
        this->msgcache_ = ss.str();
      }

      std::string msg_;
      std::size_t index_;
      std::string msgcache_;
    };

  } // namespace models
} // namespace viennashe

#endif /* VIENNASHE_MODELS_MARKOVCHAIN_EXCEPTION_HPP */

