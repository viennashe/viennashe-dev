#ifndef VIENNASHE_MATH_SPHERICAL_HARMONICS_HPP
#define VIENNASHE_MATH_SPHERICAL_HARMONICS_HPP

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


// std
#include <math.h>

// viennashe
#include "viennashe/forwards.h"

#include "viennashe/log/log.hpp"
#include "viennashe/math/log_keys.h"


/** @file viennashe/math/spherical_harmonics.hpp
    @brief Implementation of spherical harmonics plus helper functions.

    Might be replaced by functionality from ViennaMath in the future.
*/

namespace viennashe
{
  namespace math
  {


    /** @brief Compute the factorial */
    inline double factorial(long n)
    {
        double result = 1;
        for (long i=n; i>1; --i)
            result *= i;

        return result;
    }

    /** @brief Compute the double factorial, i.e. n(n-2)(n-4)... */
    inline double doublefactorial(long n)
    {
        double result = 1;
        double loopindex = n;

        while (loopindex > 1)
        {
            result *= loopindex;
            loopindex -= 2;
        }

        return result;
    }

    /** @brief Associated Legendre polynomial */
    class AssocLegendre
    {
        public:
            AssocLegendre(long n, long m) : n_(n), m_(m) {}

            double operator()(double x) const //evaluation is carried out via a recursion formula
            {
                log::info<log_AssocLegendre>() << "AssocLegendre() with n=" << n_ << ", m=" << m_ << std::endl;

                if (n_ == m_)
                    return doublefactorial(2*n_ - 1) * pow(sqrt(1.0 - x*x), n_) * ( (n_ % 2 == 0) ? 1.0 : -1.0);
                else if (m_ == (n_ - 1))
                    return x * (2.0 * n_ - 1.0) * AssocLegendre(m_, m_)(x);
                else if (n_ < m_)
                {
                    log::error() << "Error in AssocLegendre: n<m" << std::endl;
                    return 0.0;
                }
                else
                    return (  (2.0 * n_ - 1.0) * x * AssocLegendre(n_ - 1, m_)(x)
                            - (n_ + m_ - 1.0)     * AssocLegendre(n_ - 2, m_)(x) ) / (n_ - m_);
            }

        private:
            long n_;
            long m_;
    };



    /** @brief A spherical harmonic */
    class SphericalHarmonic
    {
        public:
            SphericalHarmonic(int n, int m) : n_(n), m_(m)
            {
              const double pi = 3.1415926535897932384626433832795;
              normalisation = sqrt( (2.0*n_ + 1.0) * factorial(n_ - abs(m_)) / (2.0 * pi * factorial(n_ + abs(m_))));
            }

            double operator()(double theta, double phi) const
            {
                if (m_ > 0)
                    return normalisation * AssocLegendre(n_, m_)(cos(theta)) * cos(m_ * phi);
                else if (m_ == 0)
                    return normalisation * AssocLegendre(n_, 0)(cos(theta)) / sqrt(2.0);
                else
                    return normalisation * AssocLegendre(n_, -m_)(cos(theta)) * sin(m_ * phi);
            }

        private:
          double normalisation;
          int n_;
          int m_;
    };

    /** @brief Derivative of a spherical harmonic with respect to theta */
    class SphericalHarmonic_dTheta
    {
        public:
            SphericalHarmonic_dTheta(int n, int m) : n_(n), m_(m), Y_n_m(n, m), Y_n_1_m(n-1, m) {}

            double operator()(double theta, double phi) const
            {
                if (n_ == 0)
                    return 0.0;

                if (abs(m_) < n_)
                {
                    //double sum1 = n_ * cos(theta) * SphericalHarmonic(n_, m_)(theta, phi) ;
                    //double sum2 = sqrt( (n_*n_ - m_*m_) * (2.0*n_ + 1.0) / (2.0*n_ - 1.0) ) * SphericalHarmonic(n_ - 1, m_)(theta, phi);
                    double sum1 = n_ * cos(theta) * Y_n_m(theta, phi);
                    double sum2 = sqrt( (n_*n_ - m_*m_) * (2.0*n_ + 1.0) / (2.0*n_ - 1.0) ) * Y_n_1_m(theta, phi);

                    //log::debug<log_AssocLegendre>() << sum1 << ", " << sum2 << std::endl;
                    return ( sum1 - sum2 ) / sin(theta);
                }
                else    //n == abs(m)
                    return n_ * cos(theta) * SphericalHarmonic(n_, m_)(theta, phi) / sin(theta);
            }

        private:
          int n_;
          int m_;
          SphericalHarmonic Y_n_m;    //Y_{n,m}
          SphericalHarmonic Y_n_1_m;  //Y_{n-1,m}
    };

    /** @brief Derivative of a spherical harmonic with respect to theta, alternative evaluation of theta-derivative */
    class SphericalHarmonic_dTheta_2
    {
        public:
            SphericalHarmonic_dTheta_2(int n, int m) : n_(n), m_(m) {}

            double operator()(double theta, double phi) const
            {
                if (n_ == 0)
                    return 0.0;

                const double pi = 3.1415926535897932384626433832795;
                double normalisation = sqrt( (2.0*n_ + 1.0) * factorial(n_ - abs(m_)) / (2.0 * pi * factorial(n_ + abs(m_))));

                double dP_dtheta = n_ * cos(theta) * AssocLegendre(n_, abs(m_))(cos(theta)) ;
                if (abs(m_) < n_)
                    dP_dtheta -= (n_ + abs(m_)) * AssocLegendre(n_-1, abs(m_))(cos(theta));
                dP_dtheta /= sin(theta);

                if (m_ > 0)
                    return normalisation * dP_dtheta * cos(m_ * phi);
                else if (m_ == 0)
                    return normalisation * dP_dtheta / sqrt(2.0);
                else
                    return normalisation * dP_dtheta * sin(m_ * phi);
            }

        private:
            int n_;
            int m_;
    };

    /** @brief Derivative of a spherical harmonic with respect to phi */
    class SphericalHarmonic_dPhi
    {
        public:
            SphericalHarmonic_dPhi(int n, int m) : n_(n), m_(m) {}

            double operator()(double theta, double phi) const
            {
                if (m_ == 0)
                    return 0.0;

                return m_ * SphericalHarmonic(n_, -1 * m_)(theta, phi);
            }


        private:
            int n_;
            int m_;
    };



    //Iteration over spherical harmonics indices:


    /** @brief Iteration over all spherical harmonics indices up to order L, with increasing indices l and m */
    class spherical_harmonics_iterator
    {
      public:
        spherical_harmonics_iterator(long Lmax, harmonics_iteration_type it_type)
          : L_(Lmax),
            l_( (it_type == ODD_HARMONICS_ITERATION_ID) ?  1 : 0),
            m_( (it_type == ODD_HARMONICS_ITERATION_ID) ? -1 : 0),
            iter_type_(it_type) {}

        void operator++()
        {
          if (m_ < l_)
            ++m_;
          else
          {
            switch (iter_type_)
            {
            case ALL_HARMONICS_ITERATION_ID:
              l_ += 1;
              break;
            case EVEN_HARMONICS_ITERATION_ID:
            case  ODD_HARMONICS_ITERATION_ID:
              l_ += 2;
              break;
            }
            m_ = -1 * l_;
          }
        }
        long valid() { return l_ <= L_; }
        long operator*() const { return l_*l_ + m_ + l_; }

        long index1() const { return l_; }
        long index2() const { return m_; }

      private:
        long L_;
        long l_;
        long m_;
        harmonics_iteration_type iter_type_;
      };


  } //namespace math
} //namespace viennashe
#endif
