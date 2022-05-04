#ifndef VIENNASHE_MATH_BERNOULLI_HPP
#define VIENNASHE_MATH_BERNOULLI_HPP

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


#include <math.h>

/** @file viennashe/math/bernoulli.hpp
    @brief Implementation of the Bernoulli function

    To be replaced by functionality in ViennaMath
*/

namespace viennashe
{
  namespace math
  {

    /** @brief The Bernoulli function f(x) = x / (exp(x) - 1). Avoids round-off errors near zero.*/
    template <typename NumericT>
    NumericT Bernoulli(NumericT x)
    {
      /*if (fabs(x) > 0.001)
        return x/(exp(x) - 1.0);
      else //avoids round-off errors due to exp(x) - 1
        return 1.0 / (1.0 + x / 2.0 + x*x / 6.0);*/
      if (fabs(x) > 0.01)
        return x/(exp(x) - 1.0);
      else //avoids round-off errors due to exp(x) - 1
        return 1.0 / (1.0 + x / 2.0 + x*x / 6.0 + x*x*x/24.0);
    }


    /** @brief Derivative of the Bernoulli function.  f'(x) = [exp(x) - 1 - x exp(x)] / (exp(x) - 1)^2. Avoids round-off errors near zero.*/
    template <typename NumericT>
    NumericT Bernoulli_dx(NumericT x)
    {
      if (fabs(x) > 1e-4)
      {
        double exp_x = exp(x);
        return (exp_x - 1.0 - x * exp_x) / (exp_x - 1.0) / (exp_x - 1.0);
      }
      else // replace by first-order Taylor expansion at x = 0
        return -0.5 + x / 6.0;
    }


  } //namespace math
} //namespace viennashe
#endif
