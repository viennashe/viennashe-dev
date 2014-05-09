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

#include "viennashe/math/random.hpp"
#include "viennashe/exception.hpp"

/** \file random_numbers.cpp Contains tests for the functionality provided in random.hpp (Random Number Generators)
 *  \test Tests the random number generators provided in random.hpp
 */

/** @brief Tests the std::rand based random number generator. Throws if the test fails */
inline void test_std_generator()
{
  viennashe::log::info() << "random_numbers-test: Testing std::rand based generator ..." << std::endl;

  const double values[] = {
    0.565811, 0.61093, 0.505768, 0.179647, 0.816686, 0.183472, 0.584653, 0.422156, 0.0253342,
    0.31626, 0.0612762, 0.083659, 0.976791, 0.978058, 0.873889, 0.0530768, 0.268988, 0.0923382, 0.880963,
    0.562573, 0.663514, 0.981441, 0.96191, 0.139447, 0.223442, 0.214417, 0.0119859, 0.868391, 0.331787,
    0.536112, 0.556597, 0.897598, 0.147042, 0.0623649, 0.0772447, 0.963728, 0.245837, 0.661898, 0.385884,
    0.271171, 0.978157, 0.44716, 0.35483, 0.954948, 0.425219, 0.228719, 0.00802495, 0.694207, 0.321057,
    0.888988, 0.25678, 0.984571, 0.870429, 0.218691, 0.124018, 0.0938713, 0.433108, 0.136004, 0.962262,
    0.764895, 0.672116, 0.518859, 0.662493, 0.819158, 0.581224, 0.739737, 0.782886, 0.82706, 0.401635,
    0.16877, 0.0982309, 0.379792, 0.615931, 0.453061, 0.33474, 0.0411494, 0.681779, 0.342765, 0.735357,
    0.00283636, 0.231754, 0.992137, 0.987407, 0.102183, 0.210828, 0.111425, 0.196054, 0.643936, 0.247429,
    0.158316, 0.408831, 0.919545, 0.677175, 0.0713235, 0.738702, 0.258398, 0.811061, 0.521588, 0.0854583, 0.212696
  };

  viennashe::math::std_rand_generator<double> rnd;
  rnd.seed(10);

  viennashe::log::info() << "random_numbers-test: Starting the random numbers test for 100 pseudo-random numbers. Seed = 10" << std::endl;

  for (long i = 0; i < 100; ++i)
  {
    const double r = rnd();
    if (!viennashe::testing::fuzzy_equal(r, values[i], 1e-5)) // Need fuzzy equal since we cut the numbers after the 6th post-comma pos
      throw viennashe::invalid_value_exception("random_numbers-test: 'Invalid' pseudo-random number found! ", r);
  }
} // test_std_generator

/** @brief Tests the merseinne twister 19937 RNG. Throws if the test fails*/
inline void test_merseinne_generator()
{
  viennashe::log::info() << "random_numbers-test: Testing MT 19937 based generator ..." << std::endl;

  const double values[] = {
0.07630829117,0.2273390749,0.7799187957,0.3189722276,
0.4384092249,0.9782228961,0.7234651798,0.4555849077,0.9779895162,0.308012767,0.5384958717,0.2638708407,0.5011204649,0.08674343512,0.07205113349,0.4193722107,
0.268438986,0.0159103591,0.4998824985,0.5277647912,0.6792299906,0.8688014559,0.8037390409,0.3308392481,0.380941133,0.3929423094,0.06593634444,0.6743304261,
0.2881456004,0.6723172679,0.909593527,0.6940315804,0.2133853538,0.3459729401,0.4521239561,0.9295281898,0.9312060233,0.2625837752,0.02489922685,0.7507627271,
0.6005489219,0.2548940622,0.9501294941,0.851294585,0.2303028833,0.1740527621,0.5484899166,0.7907635062,0.9091283688,0.9376293379,0.1331694445,0.4488257968,
0.5234125785,0.3846495606,0.7504098583,0.3574752405,0.6690132432,0.1933556201,0.467752859,0.1004739739,0.2048490918,0.482636906,0.4907658913,0.6137034697,
0.3723846888,0.9512483468,0.4774011539,0.3487563783,0.365890386,0.2765147036,0.8379179917,0.4465058893,0.7686475061,0.9579137282,0.3139946812,0.1795768952,
0.5726253369,0.1427157447,0.2760490526,0.08350809757,0.4528429292,0.547235124,0.3529783643,0.3102092433,0.657399467,0.2827197241,0.3703510805,0.920993543,
0.4590929756, 0.4333421632, 0.7193241196, 0.4187572289, 0.4129918288, 0.9862231447, 0.9064232644, 0.793569999, 0.1804516234, 0.3433051, 0.7411188779, 0.3003390802
  };

  viennashe::math::merseinne_twister_generator<double> rnd;
  rnd.seed(7);

  viennashe::log::info() << "random_numbers-test: Starting the random numbers test for 100 pseudo-random numbers. Seed = 7" << std::endl;

  for (long i = 0; i < 100; ++i)
  {
    const double r = rnd();
    if (!viennashe::testing::fuzzy_equal(r, values[i], 1e-7))
      throw viennashe::invalid_value_exception("random_numbers-test: 'Invalid' pseudo-random number found! ", r);
  }
} // test_merseinne_generator



/*
 * @brief A primitive test for pseudo-random numbers
 */
int main()
{
  viennashe::log::info() << "random_numbers-test: Started ..." << std::endl;

  test_std_generator();

  test_merseinne_generator();

  viennashe::log::info() << "random_numbers-test: Finished!" << std::endl;

  return (EXIT_SUCCESS);
}

