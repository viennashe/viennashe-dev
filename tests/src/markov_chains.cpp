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

#include <cstdlib>

#include "tests/src/common.hpp"

#include "viennashe/models/markovchain/reaction_rates.hpp"
#include "viennashe/models/markovchain/chain.hpp"
#include "viennashe/math/random.hpp"
#include "viennashe/models/markovchain/ssa.hpp"

/** \file markov_chains.cpp Contains tests for the markov chains, reaction rates and the SSA algorithm.
 *  \test Tests markov chains, primitive fixed reaction rates and the SSA algorithm
 */


int main()
{
  const double tol = 1e-9;

  /*
   * Expected output: "

           | to     0 | to     1 |
from     0 |        0 |       12 |
from     1 |       21 |        0 |


   States          0 |          1 |
Occupancy          1 |          0 |


dt = 0.0786663
   States          0 |          1 |
Occupancy          0 |          1 |
"
   */


  // Random number generator
  viennashe::math::std_rand_generator<double> rnd;
  rnd.seed(10);

  // Rates
  viennashe::models::const_rate r12(12);
  viennashe::models::const_rate r21(21);

  // States
  viennashe::models::state_base s1(0, "state 1");
  viennashe::models::state_base s2(1, "state 2");

  // Chain
  viennashe::models::chain mychain;
  // Allocate chain
  mychain.add_rate(s1, s2, r12);
  mychain.add_rate(s2, s1, r21);

  // Print the rate matrix
  std::cout << std::endl;
  mychain.print_rate_matrix();
  std::cout << std::endl;

  const std::size_t size = mychain.size1();
  if (size != 2 ) throw viennashe::invalid_value_exception("markov_chains-test: Size-test -> should 2, but is not ", static_cast<double>(size));

  if (std::fabs(mychain.get_rate(0,1) - r12.value()) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: rate test failed. Invalid rate: ", mychain.get_rate(0, 1));
  if (std::fabs(mychain.get_rate(1,0) - r21.value()) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: rate test failed. Invalid rate: ", mychain.get_rate(1, 0));
  if (std::fabs(mychain.get_rate(0,0)) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: rate test failed. Invalid rate: ", mychain.get_rate(0, 0));
  if (std::fabs(mychain.get_rate(1,1)) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: rate test failed. Invalid rate: ", mychain.get_rate(1, 1));


  // Solver (SSA))
  viennashe::models::ssa_solver ssa(mychain);
  // Set the equilibrium !!
  ssa.set_equilibrium(rnd);
  // Print occupancies
  std::cout << std::endl;
  mychain.print_occupancies();
  std::cout << std::endl;

  if (std::fabs(mychain.get_state(0).occupancy() - 1.0) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: The occupancy of state 0 should be 1. ", mychain.get_state(0).occupancy());
  if (std::fabs(mychain.get_state(1).occupancy()) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: The occupancy of state 1 should be 0. ", mychain.get_state(1).occupancy());

  // Solve per SSA and print occupancies
  double dt = ssa.solve(rnd);
  std::cout << std::endl << "dt = " << std::setprecision(10) << dt << std::endl;
  if (!viennashe::testing::fuzzy_equal(dt, 0.07866631598, tol))
    throw viennashe::invalid_value_exception("markov_chains-test: dt should be '0.07866631598' (tol = 1e-10) ", dt);

  mychain.print_occupancies();

  if (std::fabs(mychain.get_state(1).occupancy() - 1.0) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: The occupancy of state 1 should be 1. ", mychain.get_state(1).occupancy());
  if (std::fabs(mychain.get_state(0).occupancy()) > 1e-10)
    throw viennashe::invalid_value_exception("markov_chains-test: The occupancy of state 0 should be 0. ", mychain.get_state(0).occupancy());

  std::cout << std::endl;
  std::cout << "Test finished successfully!" << std::endl;

  return (EXIT_SUCCESS);
}

