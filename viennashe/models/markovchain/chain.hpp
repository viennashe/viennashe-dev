#ifndef VIENNASHE_MODELS_MARKOVCHAIN_CHAIN_HPP
#define VIENNASHE_MODELS_MARKOVCHAIN_CHAIN_HPP
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

// std
#include <map>

// viennashe
#include "viennashe/models/markovchain/reaction_rates.hpp"
#include "viennashe/exception.hpp"
#include "viennashe/models/markovchain/exception.hpp"
#include "viennashe/log/log.hpp"
#include "viennashe/math/linalg_util.hpp"
#include "viennashe/solvers/forwards.h"


/** @file viennashe/models/markovchain/chain.hpp
    @brief Implements the concept of a state (in a finite-state diagram) and markov chains
*/

namespace viennashe
{
  namespace models
  {

    /** @brief The basic concept of a state in Markov-Chains */
    class state_base
    {
    public:
      typedef std::size_t index_type; // The index type used to uniquely identify states

      state_base(index_type id, std::string const & name ) : id_(id), name_(name), occ_(0) { }

      state_base(state_base const & s) : id_(s.id_), name_(s.name_), occ_(s.occ_) {  }

      virtual ~state_base() { }

      /** @brief Returns the unique id of the state */
      index_type  id()   const { return id_;   }
      /** @brief Returns the name of the state for convenience and debugability */
      std::string name() const { return name_; }

      /** @brief Returns the occupancy f or p of the state. It has to be a number in [0,1] */
      double occupancy() const { return occ_; }

      /** @brief Sets the states occupancy.
       * @param o The new occupancy. Must be a number in [0,1]
       */
      void occupancy(double o)
      {
        if (o < 0.0)
          throw viennashe::invalid_value_exception("state_base.occpuancy(doule o): The occupancy must be greater or equal 0!", o);
        else if (o > 1.0)
          throw viennashe::invalid_value_exception("state_base.occpuancy(doule o): The occupancy must be less or equal 1!", o);

        occ_ = o;
      }

      /** @brief Basic interface to get a clone (uses new operator!) of the current object for storage.
       *         The caller takes ownership!
      */
      virtual state_base * clone() const { return new state_base(*this); }

    private:
      index_type  id_;
      std::string name_;
      double      occ_;
    };


    /** @brief Implementation of Markov-Chains */
    class chain
    {
    public:
      typedef state_base::index_type index_type; // The basic index type to identify states

    private:

      typedef rate_base* rate_ptr_type;

      typedef std::map<index_type, rate_ptr_type > rate_map_type;

      typedef std::map<index_type, rate_map_type > transition_map_type;

      typedef state_base*   state_ptr_type;
      typedef std::map< index_type, state_ptr_type > state_map_type;

    public:

      /** @brief Adds a single state to the chain */
      void add_state(state_base const & st)
      {
        if(states_.count(st.id()) <= 0)
        {
          chainmap_[st.id()] = rate_map_type();
          states_[st.id()]   = st.clone();
        }
      }

      /** @brief Adds a rate between to states. Additionally adds the given states if necessary */
      void add_rate(state_base const & from, state_base const & to, rate_base const & r)
      {
        this->add_rate(from.id(), to.id(), r);
        if(states_.count(from.id()) <= 0) states_[from.id()] = from.clone();
        if(states_.count(to.id()) <= 0)   states_[to.id()]   = to.clone();
      }

      /** @brief Returns a constant reference to a rate between the given states */
      rate_base const & rate(state_base const & from, state_base const & to) const
      {
        return this->rate(from.id(), to.id());
      }

      /** @brief Returns a constant reference to a rate between to states identified by their resp. ids */
      rate_base const & rate(index_type from, index_type to) const
      {
        transition_map_type::const_iterator chainmap_it = chainmap_.find(from);
        if (chainmap_it == chainmap_.end())
          throw std::runtime_error("Unknown from-state");
        rate_map_type const & frommap = chainmap_it->second;
        rate_map_type::const_iterator frommap_it = frommap.find(to);
        if (frommap_it == frommap.end())
          throw std::runtime_error("Unknown to-state");
        return *(frommap_it->second);
      }

      /** @brief Returns the actual rate (double value!) between to states */
      double get_rate(state_base const & from, state_base const & to) const
      {
        return this->get_rate(from.id(), to.id());
      }

      /** @brief Returns the row size of the rate-matrix */
      index_type size1() const { return chainmap_.size(); }
      /** @brief Returns the column size of the rate-matrix */
      index_type size2() const { return chainmap_.size(); }

      /** @brief Returns the acutal rate (double value!) between to states identified by their resp. ids */
      double get_rate(index_type from, index_type to) const
      {
        transition_map_type::const_iterator chainmap_it = chainmap_.find(from);
        if (chainmap_it == chainmap_.end())
          throw std::runtime_error("Unknown from-state");
        rate_map_type const & frommap = chainmap_it->second;
        rate_map_type::const_iterator frommap_it = frommap.find(to);
        if (frommap_it == frommap.end())
          return 0;
        return frommap_it->second->value();
      }

      /** @brief Prints the rate matrix onto screen (log::info). Usefull for debugging */
      void print_rate_matrix() const
      {
        std::vector<index_type> ids(chainmap_.size());
        std::size_t i = 0;
        viennashe::log::info() << std::setw(10) << " " <<  " | ";
        for (transition_map_type::const_iterator cit = chainmap_.begin();
            cit != chainmap_.end(); ++cit)
        {
          viennashe::log::info() << "to " << std::setw(5) << cit->first << " | ";
          ids[i++] = cit->first;
        }
        viennashe::log::info() << std::endl;

        for (transition_map_type::const_iterator cit = chainmap_.begin();
            cit != chainmap_.end(); ++cit)
        {
          viennashe::log::info() << "from " << std::setw(5) << cit->first << " | ";
          for (index_type j = 0; j < ids.size(); ++j)
          {
            bool found = false;
            for (rate_map_type::const_iterator rit = cit->second.begin();
                 rit != cit->second.end(); ++rit)
            {
              if (rit->first == j)
              {
                viennashe::log::info() << std::setw(8) << rit->second->value() << " | ";
                found = true;
                break;
              }
            }
            if (! found)
            {
              viennashe::log::info() << std::setw(8) << 0 << " | ";
            }
          }
          viennashe::log::info() << std::endl;
        }
      } // print_rate_matrix()

      /**
       * @brief Fills the given rate_matrix
       * @param rate_matrix An ublas compatible matrix. The matrix needs to be properly allocated!
       */
      template < typename MatrixT >
      void get_rate_matrix(MatrixT & rate_matrix) const
      {
        if (rate_matrix.size1() != this->size1())
          throw viennashe::models::invalid_parameter_exception("chain.get_rate_matrix(): Invalid matrix row size!");
        if (rate_matrix.size2() != this->size2())
          throw viennashe::models::invalid_parameter_exception("chain.get_rate_matrix(): Invalid matrix column size!");

        for (transition_map_type::const_iterator cit = chainmap_.begin();
            cit != chainmap_.end(); ++cit)
        {
          for (rate_map_type::const_iterator rit = cit->second.begin();
               rit != cit->second.end(); ++rit)
          {
            rate_matrix(cit->first, rit->first) = rit->second->value();
          }
        }
      } // get_rate_matrix()

      /** @brief Prints the occupancies of each state in the chain onto screen (log::info). Usefull for debugging */
      void print_occupancies() const
      {
        viennashe::log::info() << std::setw(10) << "States ";
        for (state_map_type::const_iterator sit = states_.begin();
             sit != states_.end(); ++sit)
        {
           viennashe::log::info() << std::setw(10) << sit->first << " | ";
        }
        viennashe::log::info() << std::endl;
        viennashe::log::info() << "Occupancy ";
        for (state_map_type::const_iterator sit = states_.begin();
             sit != states_.end(); ++sit)
        {
           viennashe::log::info() << std::setw(10) << sit->second->occupancy() << " | ";
        }
        viennashe::log::info() << std::endl;

      }

      /** @brief Returns true if a state identified by its id is in the chain */
      bool has_state(index_type idx) const { return (states_.count(idx) > 0 ); }

      state_base       & get_state(index_type idx)
      {
        if (states_.count(idx) <= 0)
          throw std::runtime_error("State not registered!");
        return *(states_[idx]);
      }
      state_base const & get_state(index_type idx) const
      {
        state_map_type::const_iterator states_it = states_.find(idx);
        if (states_it == states_.end())
          throw std::runtime_error("Unknown state");
        return *(states_it->second);
      }

      ~chain()
      {
        // clean up pointers in transition map:
        for (transition_map_type::const_iterator cit = chainmap_.begin(); cit != chainmap_.end(); ++cit)
          for (rate_map_type::const_iterator rit = cit->second.begin(); rit != cit->second.end(); ++rit)
            delete rit->second;

        // clean up pointers in state map:
        for (state_map_type::const_iterator sit = states_.begin(); sit != states_.end(); ++sit)
            delete sit->second;
      }

    protected:

      void add_rate(index_type from, index_type to, rate_base const & r)
      {
        if (chainmap_.count(from) > 0)
        {
          rate_map_type & frommap = chainmap_[from];
          frommap[to] = r.clone();
        }
        else
        {
          rate_map_type frommap;
          chainmap_[from] = frommap;
          delete chainmap_[from][to]; //safe to call
          chainmap_[from][to] = r.clone();
        }
      }

    private:
      transition_map_type chainmap_;
      state_map_type      states_;
    };

    /**
     * @brief Solves the time averaged set of equilibrium (dt->infty) equations for the given Markov-Chain and sets the respective occupancies
     *        The occupancies f_i will be values between zero and one, where sum_i f_i = 1
     * @param c A valid reference to the chain. The internal state of the given object will be changed
     */
    inline void set_chain_to_equilibrium_expectation(chain & c)
    {
      typedef viennashe::math::dense_matrix<double> MatrixType;
      typedef std::vector<double>  VectorType;

      typedef MatrixType::size_type size_type;

      size_type num_rows = c.size1();

      MatrixType A(num_rows, num_rows);
      VectorType b(num_rows);
      VectorType x(num_rows);

      // Set RHS to zero
      for (size_type i = 0; i < num_rows; ++i)
      {
        b[i] = 0.0;
        for (size_type j = 0; j < num_rows; ++j)
          A(i,j) = 0.0;
      }

      // Fill dense matrix
      for (size_type i = 0; i < num_rows; ++i)
      {
        if (! c.has_state(i))
          throw viennashe::models::invalid_state_exception("set_chain_to_equilibrium_expectation(): State does not exist!", i);

        for (size_type j = 0; j < num_rows; ++j)
        {
          const double kij = c.get_rate(i, j);
          A(i,i) -= kij;
          A(j,i) += kij;
        }
      }
      // Reset last row to contain the normalization condition sum_i f_i = 1
      b[num_rows-1] = 1.0;
      for (size_type j = 0; j < num_rows; ++j)
        A(num_rows-1, j) = 1.0;

      // DEBUG
      /*
      for (size_type i = 0; i < num_rows; ++i)
      {
        for (size_type j = 0; j < num_rows; ++j)
        {
          viennashe::log::debug() << std::setw(7) << A(i,j) << " ";
        }
        viennashe::log::debug() << std::endl;
      }
      */

      // Solve
      x = viennashe::solvers::solve(A, b);
      // Check solution
      double total = 0;
      for (size_type i = 0; i < num_rows; ++i)
      {
        total += x[i];
      }
      // FUZZY CHECK
      if (total - 1e-4 > 1.0 || total + 1e-4 < 0.0)
        throw viennashe::models::model_evaluation_exception("set_chain_to_equilibrium_expectation(): The sum of all occupancies is not 1.0!");
      // Set solution
      for (size_type i = 0; i < num_rows; ++i)
      {
        c.get_state(i).occupancy(x[i]);
      }

    } // set_chain_to_equilibrium_expectation

    /**
     * @brief Solves the time averaged set of equations (including time derivatives) for the occupancies f_i in a
     *        given Markov-Chain and sets the resp. occupancies.
     *        The occupancies f_i will be values between zero and one, where sum_i f_i = 1.
     *        In order to solve the time dependent equations the backward-euler method is used
     * @param c A valid reference to the chain. The internal state of the given object will be changed
     * @param dt The finite time step (in seconds)
     */
    inline void solve_for_expectation_occupancies(chain & c, double dt)
    {
      if(dt <= 0.0) throw viennashe::invalid_value_exception("solve_for_expectation_occupancies(): dt must be greater zero!", dt);

      typedef viennashe::math::dense_matrix<double> MatrixType;
      typedef std::vector<double>  VectorType;

      typedef MatrixType::size_type size_type;

      const size_type num_rows = c.size1();
      const double rtime = 1.0/dt;

      MatrixType A(num_rows, num_rows);
      VectorType b(num_rows);
      VectorType x(num_rows);

      // Set RHS to f_old/dt
      for (size_type i = 0; i < num_rows; ++i)
      {
        if (! c.has_state(i))
          throw viennashe::models::invalid_state_exception("solve_for_expectation_occupancies(): State does not exist!", i);

        b[i] = rtime * c.get_state(i).occupancy(); // f_old/dt
        for (size_type j = 0; j < num_rows; ++j)
          A(i,j) = 0.0;
      }

      // Fill dense matrix
      for (size_type i = 0; i < num_rows; ++i)
      {
        if (! c.has_state(i))
          throw viennashe::models::invalid_state_exception("solve_for_expectation_occupancies(): State does not exist!", i);

        A(i,i) += rtime;

        for (size_type j = 0; j < num_rows; ++j)
        {
          const double kij = c.get_rate(i, j);
          A(i,i) -= kij;
          A(j,i) += kij;
        }
      }
      // Reset last row to contain the normalization condition sum_i f_i = 1
      b[num_rows-1] = 1.0;
      for (size_type j = 0; j < num_rows; ++j)
        A(num_rows-1, j) = 1.0;

      // DEBUG
      /*
      for (size_type i = 0; i < num_rows; ++i)
      {
        for (size_type j = 0; j < num_rows; ++j)
        {
          viennashe::log::debug() << std::setw(7) << A(i,j) << " ";
        }
        viennashe::log::debug() << std::endl;
      }
      */

      // Solve
      x = viennashe::solvers::solve(A, b);
      // Check solution
      double total = 0;
      for (size_type i = 0; i < num_rows; ++i)
      {
        total += x[i];
      }
      // FUZZY CHECK
      if (total - 1e-4 > 1.0 || total + 1e-4 < 0.0)
        throw viennashe::models::model_evaluation_exception("solve_for_expectation_occupancies(): The sum of all occupancies is not 1.0!");
      // Set solution
      for (size_type i = 0; i < num_rows; ++i)
      {
        c.get_state(i).occupancy(x[i]);
      }

    } // solve_for_expectation_occupancies()

  } // namespace models
} // namespace viennashe


#endif /* VIENNASHE_MODELS_MARKOVCHAIN_CHAIN_HPP */

