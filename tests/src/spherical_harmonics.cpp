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

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>

#include "tests/src/common.hpp"

#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/she/harmonics_coupling.hpp"

#include <math.h>


/** \file spherical_harmonics.cpp Contains test for the spherical harmonics implementation
 *  \test Tests the orthogonality of spherical harmonics
 */


//up to which order assoc. Legendre functions should be checked:
#define LEGENDRE_DEGREE  7

//up to which order spherical harmonics should be checked (values of 5 and above may take very long...)
#define SH_DEGREE  3

//up to which order coupling matrices shall be computed (values of 5 and above may take very long...)
#define COUPLING_DEGREE  3

class AssocLegendreTester
{
  public:
    AssocLegendreTester(long l, long m,
                        long lprime, long mprime) : A(l, m), B(lprime, mprime) {}

    double operator()(double x) const
    {
      return A(x) * B(x);
    }


  private:
    viennashe::math::AssocLegendre A;
    viennashe::math::AssocLegendre B;
};


class SphericalHarmonicTester
{
  public:
    SphericalHarmonicTester(int l, int m,
                            int lprime, int mprime) : A(l, m), B(lprime, mprime) {}

    double operator()(double theta, double phi) const
    {
      return A(theta, phi) * B(theta, phi) * sin(theta);
    }


  private:
    viennashe::math::SphericalHarmonic A;
    viennashe::math::SphericalHarmonic B;
};


inline void fuzzy_check(double is, double should, std::string message)
{
  const double tol = 1e-7;  //tolerance
  if (!viennashe::testing::fuzzy_equal(is, should, tol)
      && std::max(std::abs(is), std::abs(should)) > 1e-4)
  {
    std::cerr << "Failed at: " << message << std::endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    //std::cout << "Success!" << std::endl;
  }
}


int main()
{
  const double pi = 3.1415926535897932384626433832795;
  viennashe::math::IntAdaptiveRule< viennashe::math::IntGauss<5> > integration_rule;
  viennashe::math::SphericalHarmonic Y_00(0, 0);

  std::cout << "Test 1: Value of Y_00" << std::endl;
  fuzzy_check(Y_00(0, 0), 1.0 / sqrt(4.0 * pi), "Value of Y_00" );

  std::cout << "Test 2: Orthonormality of associated Legendre functions up to degree " << LEGENDRE_DEGREE << std::endl;
  for (long l=0; l <= LEGENDRE_DEGREE; ++l)
  {
    for (long m=0; m<=l; ++m)
    {
      //checknig fixed second index only
      for (long lprime=0; lprime <= LEGENDRE_DEGREE; ++lprime)
      {
        if (m > lprime)
          continue;

        std::stringstream ss;
        ss << "Testing P(" << l << "," << m << ") "
            << "with P(" << lprime << "," << m << ")" << std::endl;

        double expected_result = 0.0;
        if (l == lprime)
          expected_result = 2.0 * viennashe::math::factorial(l+m)
                          / ( (2.0*l + 1.0) * viennashe::math::factorial(l-m) );

        double is_result = viennashe::math::integrate(AssocLegendreTester(l, m, lprime, m),
                                                      -1.0, 1.0,
                                                      integration_rule);
        fuzzy_check(is_result, expected_result, ss.str());
      } //for lprime

    } //for m
  } // for l

  std::cout << "Test 3: Orthonormality of Spherical Harmonics up to degree " << SH_DEGREE << std::endl;
  for (int l=0; l <= SH_DEGREE; ++l)
  {
    for (int m=-l; m<=l; ++m)
    {
      for (int lprime=0; lprime <= SH_DEGREE; ++lprime)
      {
        for (int mprime=-lprime; mprime <= lprime; ++mprime)
        {
          std::stringstream ss;
          ss << "Testing Y(" << l << "," << m << ") "
             << "with Y(" << lprime << "," << mprime << ")" << std::endl;

          double expected_result = 0.0;
          if ( (l == lprime) && (m == mprime) )
            expected_result = 1.0;

          double is_result = viennashe::math::integrate2D(0.0, pi,
                                                          0.0, 2.0 * pi,
                                                          SphericalHarmonicTester(l, m, lprime, mprime),
                                                          integration_rule);
          fuzzy_check(is_result, expected_result, ss.str());
        } //for mprime
      } //for lprime
    } //for m
  } // for l



  std::cout << "Test 4: Coupling matrices up to degree " << COUPLING_DEGREE << " (this might take a while...)" << std::endl;
  for (int l=0; l <= COUPLING_DEGREE; ++l)
  {
    for (int m=-l; m<=l; ++m)
    {
      for (int lprime=0; lprime <= COUPLING_DEGREE; ++lprime)
      {
        for (int mprime=-lprime; mprime <= lprime; ++mprime)
        {
          viennashe::she::vIntegrand_x integrand_a_x(l, m, lprime, mprime);
          viennashe::she::vIntegrand_y integrand_a_y(l, m, lprime, mprime);
          viennashe::she::vIntegrand_z integrand_a_z(l, m, lprime, mprime);

          viennashe::she::GammaIntegrand_x integrand_b_x(l, m, lprime, mprime);
          viennashe::she::GammaIntegrand_y integrand_b_y(l, m, lprime, mprime);
          viennashe::she::GammaIntegrand_z integrand_b_z(l, m, lprime, mprime);

          double result_a_x = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_a_x, integration_rule);
          double result_a_y = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_a_y, integration_rule);
          double result_a_z = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_a_z, integration_rule);

          double result_b_x = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_b_x, integration_rule);
          double result_b_y = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_b_y, integration_rule);
          double result_b_z = viennashe::math::integrate2D(0.0, pi, 0.0, 2.0 * pi, integrand_b_z, integration_rule);

          if (std::abs(l - lprime) > 1)
          {
            fuzzy_check(result_a_x, 0, "result_a_x");
            fuzzy_check(result_a_y, 0, "result_a_y");
            fuzzy_check(result_a_z, 0, "result_a_z");

            fuzzy_check(result_b_x, 0, "result_b_x");
            fuzzy_check(result_b_y, 0, "result_b_y");
            fuzzy_check(result_b_z, 0, "result_b_z");
          }

        } //for mprime
      } //for lprime
    } //for m
  } // for l


  std::cout << "*******************************" << std::endl;
  std::cout << "* Test finished successfully! *" << std::endl;
  std::cout << "*******************************" << std::endl;

} //main()

