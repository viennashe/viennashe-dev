#ifndef VIENNASHE_MODELS_MARKOVCHAIN_SSA_HPP
#define VIENNASHE_MODELS_MARKOVCHAIN_SSA_HPP
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

#include "viennashe/forwards.h"

#include "viennashe/models/markovchain/exception.hpp"
#include "viennashe/models/exception.hpp"

#include "viennashe/models/markovchain/reaction_rates.hpp"
#include "viennashe/models/markovchain/chain.hpp"
#include "viennashe/math/random.hpp"

#include "viennashe/math/linalg_util.hpp"
#include "viennashe/solvers/forwards.h"

/** @file viennashe/models/markovchain/ssa.hpp
    @brief Contains an impelemtation of the SSA algorithm for Markov-Chains (cf. viennashe/models/markovchain/chain.hpp)
 */

namespace viennashe
{
  namespace models
  {

    /** @brief The SSA (Monte Carlo) solver for Markov-Chains
     *         The given chain is being referenced
     *         and the occupancies p_i of the states are being changed.
     *         The resp. occupancies can either be 1.0 (occupied) or 0 (unoccupied)!
     *         The solver additionally guarantees that sum_i p_i = 1.0, i.e. only one state is occupied at a time !
     */
    class ssa_solver
    {
    public:
      typedef viennashe::models::chain::index_type index_type;

      /** @brief CTOR. */
      ssa_solver(viennashe::models::chain & c) : chain_(c), dt_(0) { }

      /** @brief Returns the last time step dt in seconds */
      double delta_t() const { return dt_; }

      /**
       * @brief Sets the equilibrium occupancies
       * @param rnd A (pseudo-) random number generator
       */
      void set_equilibrium(viennashe::math::rand_generator_base<double> & rnd)
      {
        // This should fill all the occupancies with values between 0 and 1 where sum_i f_i = 1 !!
        viennashe::models::set_chain_to_equilibrium_expectation(chain_);
        // Debug
        //chain_.print_occupancies();

        const double random_value = rnd(); // [0,1) ...

        double sum = 0;
        index_type i = 0;
        for(i = 0; i < chain_.size1(); ++i)
        {
          sum += chain_.get_state(i).occupancy(); // Now this should be the probability 0 <= p_i <= 1 !!
          if (random_value <= sum)
            break;
        }
        // Set all occupancies to zero
        for(index_type j = 0; j < chain_.size1(); ++j)
          chain_.get_state(j).occupancy(0.0);
        // Set the lucky occupancy to 1
        chain_.get_state(i).occupancy(1.0);
      }

      /**
       * @brief Determines the time until the next transition event and advances the system to this state
       * @param rnd A (pseudo-) random number generator
       * @return The time (delta t) until the next event. The internal state of the chain is being changed.
       */
      double solve(viennashe::math::rand_generator_base<double> & rnd)
      {
        double dt = 0;
        const double a0 = this->sum_rates();

        if (a0 <= 0.0)
          throw viennashe::models::model_evaluation_exception("ssa.solve(): The sum of rates is zero or negative!");

        double r1 = 1.0 - rnd();
        { while (!r1) r1 = 1.0 - rnd(); } // Exclude 0.0
        const double r2 = rnd();

        dt = 1.0/a0 * std::log( 1.0/r1 );

        double sum = 0;
        index_type i = 0;
        index_type j = 0;
        for (i = 0; i < chain_.size1(); ++i)
        {
          if(chain_.get_state(i).occupancy() != 1.0) continue;

          for (j = 0; j < chain_.size2(); ++j)
          {
            const double kij = chain_.get_rate(i, j);

            if (kij <= 0.0) continue; // skip non-existent rates

            const double pij = kij / a0;
            sum += pij;
            if (r2 <= sum) break; // we found the transition
          }
          if (r2 <= sum) break; // we found the transition
        }

        if (i >= chain_.size1())
          throw viennashe::models::model_evaluation_exception("ssa.solve(): The SSA algorithm failed to find the next state!");

        chain_.get_state(j).occupancy(1.0);
        chain_.get_state(i).occupancy(0.0);

        this->dt_ = dt; // Cache delta time
        return this->delta_t();
      }

      viennashe::models::chain const &  chain() const { return this->chain_;  }

    private:

      double sum_rates() const
      {
        double sum = 0;
        for (index_type i = 0; i < chain_.size1(); ++i)
        {
          if(chain_.get_state(i).occupancy() != 1.0) continue;

          for (index_type j = 0; j < chain_.size2(); ++j)
          {
            const double kij = chain_.get_rate(i, j);
            if (kij <= 0.0) continue;
            sum += kij;
          }
        }
        return sum;
      }

      viennashe::models::chain & chain_;
      double dt_;
    };



  } // namespace models
} // namespace viennashe


#endif /* VIENNASHE_MODELS_MARKOVCHAIN_SSA_HPP */

