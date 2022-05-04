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



#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>

#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/she/harmonics_coupling.hpp"

#include <math.h>

/** \file spherical_harmonics_iter.cpp Tests the spherical harmonics iterator
 *  \test Tests the spherical harmonics iterator
 */


/**
 * @brief Tests the SH iterators. Returns false if the tests fail
 * @param iter The iterator
 * @param ref_values Vector of reference values
 * @param ref_valid Vector containing bools, if true the reference value is valid
 */
inline bool test_iterator(viennashe::math::spherical_harmonics_iterator iter,
                          std::vector<long> ref_values,
                          std::vector<bool> ref_valid)
{
  for (std::size_t i=0; i<ref_values.size(); ++i)
  {
    if (ref_values[i] != *iter)
    {
      std::cerr << "Error after " << i << "increments: Expected value " << ref_values[i] << ", but got " << *iter << std::endl;
      return false;
    }
    if (ref_valid[i] != iter.valid())
    {
      std::cerr << "Error after " << i << "increments: Expected valid status " << ref_valid[i] << ", but got " << iter.valid() << std::endl;
      return false;
    }
    ++iter;
  }
  return true;
}

int main()
{
  viennashe::math::spherical_harmonics_iterator iter_all(3, viennashe::math::ALL_HARMONICS_ITERATION_ID);
  viennashe::math::spherical_harmonics_iterator iter_even(3, viennashe::math::EVEN_HARMONICS_ITERATION_ID);
  viennashe::math::spherical_harmonics_iterator iter_odd(3, viennashe::math::ODD_HARMONICS_ITERATION_ID);

  bool ok = false;

  std::cout << "--- Iteration over all harmonics: ---" << std::endl;
  std::vector<bool> ref_all_valid(17, true);
  ref_all_valid[16] = false;
  std::vector<long> ref_all_values(17);
  for (std::size_t i=0; i<17; ++i)
    ref_all_values[i] = long(i);

  ok = test_iterator(iter_all, ref_all_values, ref_all_valid);
  if (ok == false) return (EXIT_FAILURE);

  std::cout << "--- Iteration over even harmonics: ---" << std::endl;
  std::vector<bool> ref_even_valid(7, true);
  ref_even_valid[6] = false;
  std::vector<long> ref_even_values(7);
  ref_even_values[1] = 4;
  ref_even_values[2] = 5;
  ref_even_values[3] = 6;
  ref_even_values[4] = 7;
  ref_even_values[5] = 8;
  ref_even_values[6] = 16;
  ok = test_iterator(iter_even, ref_even_values, ref_even_valid);
  if (ok == false) return (EXIT_FAILURE);

  std::cout << "--- Iteration over odd harmonics: ---" << std::endl;
  std::vector<bool> ref_odd_valid(11, true);
  ref_odd_valid[10] = false;
  std::vector<long> ref_odd_values(11);
  ref_odd_values[0] = 1;
  ref_odd_values[1] = 2;
  ref_odd_values[2] = 3;
  ref_odd_values[3] = 9;
  ref_odd_values[4] = 10;
  ref_odd_values[5] = 11;
  ref_odd_values[6] = 12;
  ref_odd_values[7] = 13;
  ref_odd_values[8] = 14;
  ref_odd_values[9] = 15;
  ref_odd_values[10] = 25;
  ok = test_iterator(iter_odd, ref_odd_values, ref_odd_valid);
  if (ok == false) return (EXIT_FAILURE);

  std::cout << "*******************************" << std::endl;
  std::cout << "* Test finished successfully! *" << std::endl;
  std::cout << "*******************************" << std::endl;

  return (EXIT_SUCCESS);

} //main()

