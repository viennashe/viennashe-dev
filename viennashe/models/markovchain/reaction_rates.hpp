#ifndef VIENNASHE_MODELS_MARKOVCHAIN_REACTION_RATES_HPP
#define VIENNASHE_MODELS_MARKOVCHAIN_REACTION_RATES_HPP
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

/** @file viennashe/models/markovchain/reaction_rates.hpp
    @brief Contains the basic reaction rate interface.
 */

namespace viennashe
{
  namespace models
  {

    /** @brief The basic rate interface */
    struct rate_base
    {
      typedef double value_type;

      virtual ~rate_base() { }

      /**
       * The main functor interface
       * @return The rate in 1/s
       */
      virtual value_type value() const = 0;

      /**
       * Generates a clone of itself (deep copy; uses new). The caller takes ownership
       * @return A pointer to a deep copy of itself
       */
      virtual rate_base * clone() const = 0;

      /**
       * A wrapper to value()
       * @return The same as the virtual abstract method value()
       */
      virtual value_type operator()() const { return this->value(); }

    };

    /** @brief A simple constant rate */
    struct const_rate : public rate_base
    {
      typedef rate_base::value_type value_type;

      const_rate(value_type rate) : rate_(rate) { }
      ~const_rate() { }

      virtual value_type value() const { return rate_; }

      virtual rate_base * clone() const
      {
        return new const_rate(this->rate_);
      }

    private:
      value_type rate_;
    };


  } // namespace models
} // namespace viennashe



#endif /* VIENNASHE_MODELS_MARKOVCHAIN_REACTION_RATES_HPP */

