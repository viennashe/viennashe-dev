#ifndef VIENNASHE_MATH_INTEGRATOR_HPP
#define VIENNASHE_MATH_INTEGRATOR_HPP

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

///// PURPOSE: Provide integration routines in 1D and 2D
//// Nomenclature: \int_a^b f(x) dx
//// Attention: This header file should be used within SHE2D only. A (much) better stand-alone package can be found in another svn-repository

#include <math.h>
#include <iostream>

#include "viennashe/forwards.h"

#include "viennashe/log/log.hpp"
#include "viennashe/math/log_keys.h"

/** @file viennashe/math/integrator.hpp
    @brief Implementation of numerical integration routines.

    Will be replaced by functionality from ViennaMath in the future.
*/

namespace viennashe
{
  /** @brief Namespace containing mathematical tools */
  namespace math
  {


    //Tags:
    /** @brief Tag for the midpoint rule */
    class IntMidPointRule {};
    /** @brief Tag for the trapezoidal rule */
    class IntTrapezoidalRule {};
    /** @brief Tag for the simpson rule */
    class IntSimpsonRule {};

    /** @brief Tag for an adaptive integration rule
    *
    * The strategy is to subdivide the integration interval into two intervals of equal length and perform adaptive integration on both until a certain accuracy is reached.
    *
    * @tparam IntRule Any non-adaptive integration rule
    */
    template <typename IntRule>
    class IntAdaptiveRule {};

    /** @brief Gaussian integration rule
    *
    * Rule is accurate for polynomials up to degree 2 * int_points - 1
    *
    * @param int_points Number of integration points
    */
    template <int int_points>
    class IntGauss {};


    /** @brief Binds the first argument of a function of functor to a certain value
    *
    * @tparam T type of the object for which the first argument is to be bound
    */
    template <typename T>
    class BindFirst
    {
        public:
            BindFirst(T const & func, double value) : func_(func), value_(value) {}

            double operator()(double x) const { return func_(value_, x); }
            double operator()(double x1, double x2) const { return func_(value_, x1, x2); }
            double operator()(double x1, double x2, double x3) const { return func_(value_, x1, x2, x3); }
            double operator()(double x1, double x2, double x3, double x4) const { return func_(value_, x1, x2, x3, x4); }

            operator double() const { return func_(value_); }

        private:
            T const & func_;
            double value_;
    };

    /** @brief Binds the second argument of a function or functor to a certain value
    *
    * @tparam T type of the object for which the second argument is to be bound
    */
    template <typename T>
    class BindSecond
    {
        public:
            BindSecond(T const & func, double value) : func_(func), value_(value) {}

            double operator()(double x) const { return func_(x, value_); }
            double operator()(double x1, double x2) const { return func_(x1, value_, x2); }
            double operator()(double x1, double x2, double x3) const { return func_(x1, value_, x2, x3); }
            double operator()(double x1, double x2, double x3, double x4) const { return func_(x1, value_, x2, x3, x4); }

            //operator double() const { return func_(value_); }     //meaningless in this case!

        private:
            T const & func_;
            double value_;
    };



    /** @brief Integration of a function using the mid point rule */
    template <typename T>
    double integrate(T const & f, double a, double b, IntMidPointRule)
    {
        //log::debug<log_integrate>() << "Midpoint: " << f( (a+b)/2.0 ) << ", " << a << "," << b << std::endl;
        return ( f( (a+b)/2.0 ) ) * fabs(b-a);
    }

    /** @brief Integration of a function using the trapezoidal rule */
    template <typename T>
    double integrate(T const & f, double a, double b, IntTrapezoidalRule)
    {
        return ( f(a) + f(b) ) * fabs(b-a) / 2.0;
    }

    /** @brief Integration of a function using the simpson rule */
    template <typename T>
    double integrate(T const & f, double a, double b, IntSimpsonRule)
    {
        return (f( a ) + 4.0 * f( (a+b)/2.0 ) + f(b)) * fabs(b-a) / 6.0;
    }

    /** @brief Integration of a function using the Gauss rule, 1 evaluation node */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<1>)
    {
        return f( (a + b) / 2.0 ) * fabs(b-a);
    }

    /** @brief Integration of a function using the Gauss rule, 2 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<2>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) * 0.288675135;
        //log::debug<log_integrate>()  << f( xmid - xdiff ) << ", " << f( xmid + xdiff ) << std::endl;
        return (f( xmid - xdiff ) + f( xmid + xdiff )) * fabs(b-a) / 2.0;
    }

    /** @brief Integration of a function using the Gauss rule, 3 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<3>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) * 0.387298335;
        return ( f(xmid - xdiff ) * 5.0/9.0
              + f(xmid) * 8.0/9.0
              + f(xmid + xdiff ) * 5.0/9.0) * fabs(b-a) / 2.0;
    }

    /** @brief Integration of a function using the Gauss rule, 4 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<4>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.8611363115940525752239465 ) * 0.3478548451374538573730639
              + f(xmid - xdiff * 0.3399810435848562648026658 ) * 0.6521451548625461426269361
              + f(xmid + xdiff * 0.3399810435848562648026658 ) * 0.6521451548625461426269361
              + f(xmid + xdiff * 0.8611363115940525752239465 ) * 0.3478548451374538573730639) * xdiff;
    }

    /** @brief Integration of a function using the Gauss rule, 5 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<5>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.90617985 ) * 0.23692689
              + f(xmid - xdiff * 0.53846931 ) * 0.47862867
              + f(xmid                      ) * 0.56888888888889
              + f(xmid + xdiff * 0.53846931 ) * 0.47862867
              + f(xmid + xdiff * 0.90617985 ) * 0.23692689) * xdiff;
    }


    /** @brief Integration of a function using the Gauss rule, 6 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<6>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.9324695142031520278123016 ) * 0.1713244923791703450402961
              + f(xmid - xdiff * 0.6612093864662645136613996 ) * 0.3607615730481386075698335
              + f(xmid - xdiff * 0.2386191860831969086305017 ) * 0.4679139345726910473898703
              + f(xmid + xdiff * 0.2386191860831969086305017 ) * 0.4679139345726910473898703
              + f(xmid + xdiff * 0.6612093864662645136613996 ) * 0.3607615730481386075698335
              + f(xmid + xdiff * 0.9324695142031520278123016 ) * 0.1713244923791703450402961) * xdiff;
    }

    /** @brief Integration of a function using the Gauss rule, 7 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<7>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.9491079123427585245261897 ) * 0.1294849661688696932706114
              + f(xmid - xdiff * 0.7415311855993944398638648 ) * 0.2797053914892766679014678
              + f(xmid - xdiff * 0.4058451513773971669066064 ) * 0.3818300505051189449503698
              + f(xmid                                       ) * 0.4179591836734693877551020
              + f(xmid + xdiff * 0.4058451513773971669066064 ) * 0.3818300505051189449503698
              + f(xmid + xdiff * 0.7415311855993944398638648 ) * 0.2797053914892766679014678
              + f(xmid + xdiff * 0.9491079123427585245261897 ) * 0.1294849661688696932706114) * xdiff;
    }


    /** @brief Integration of a function using the Gauss rule, 8 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<8>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.96028986 ) * 0.10122854
              + f(xmid - xdiff * 0.79666648 ) * 0.22238103
              + f(xmid - xdiff * 0.52553241 ) * 0.31370665
              + f(xmid - xdiff * 0.18343464 ) * 0.36268378
              + f(xmid + xdiff * 0.18343464 ) * 0.36268378
              + f(xmid + xdiff * 0.52553241 ) * 0.31370665
              + f(xmid + xdiff * 0.79666648 ) * 0.22238103
              + f(xmid + xdiff * 0.96028986 ) * 0.10122854) * xdiff;
    }

    /** @brief Integration of a function using the Gauss rule, 9 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<9>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.9681602395076260898355762 ) * 0.0812743883615744119718922
              + f(xmid - xdiff * 0.8360311073266357942994298 ) * 0.1806481606948574040584720
              + f(xmid - xdiff * 0.6133714327005903973087020 ) * 0.2606106964029354623187429
              + f(xmid - xdiff * 0.3242534234038089290385380 ) * 0.3123470770400028400686304
              + f(xmid                                       ) * 0.3302393550012597631645251
              + f(xmid + xdiff * 0.3242534234038089290385380 ) * 0.3123470770400028400686304
              + f(xmid + xdiff * 0.6133714327005903973087020 ) * 0.2606106964029354623187429
              + f(xmid + xdiff * 0.8360311073266357942994298 ) * 0.1806481606948574040584720
              + f(xmid + xdiff * 0.9681602395076260898355762 ) * 0.0812743883615744119718922) * xdiff;
    }

    /** @brief Integration of a function using the Gauss rule, 10 evaluation nodes */
    template <typename T>
    double integrate(T const & f, double a, double b, IntGauss<10>)
    {
        double xmid  = (a + b) / 2.0;
        double xdiff = (b - a) / 2.0;

        return ( f(xmid - xdiff * 0.9739065285171717200779640 ) * 0.0666713443086881375935688
              + f(xmid - xdiff * 0.8650633666889845107320967 ) * 0.1494513491505805931457763
              + f(xmid - xdiff * 0.6794095682990244062343274 ) * 0.2190863625159820439955349
              + f(xmid - xdiff * 0.4333953941292471907992659 ) * 0.2692667193099963550912269
              + f(xmid - xdiff * 0.1488743389816312108848260 ) * 0.2955242247147528701738930
              + f(xmid + xdiff * 0.1488743389816312108848260 ) * 0.2955242247147528701738930
              + f(xmid + xdiff * 0.4333953941292471907992659 ) * 0.2692667193099963550912269
              + f(xmid + xdiff * 0.6794095682990244062343274 ) * 0.2190863625159820439955349
              + f(xmid + xdiff * 0.8650633666889845107320967 ) * 0.1494513491505805931457763
              + f(xmid + xdiff * 0.9739065285171717200779640 ) * 0.0666713443086881375935688) * xdiff;
    }


    /** @brief Integration of a function using an adaptive integration rule based on bisection of the integration interval */
    template <typename T, typename U>
    double integrate(T const & f, double a, double b, IntAdaptiveRule<U> adaptiverule, long recursion_depth = 0)
    {
        double direct = integrate(f, a, b, U());

        //log::debug<log_integrate>()  << "Calling adaptive with a=" << a << ", b=" << b << " and recursion_depth=" << recursion_depth << std::endl;

        if (recursion_depth == 10)
        {
            //double half_left = integrate(f, a, (a+b)/2.0, U());
            //double half_right = integrate(f, (a+b)/2.0, b, U());

            //log::warn<log_integrate>()  << "Recursion depth reached. a: " << a << ", b: " << b << ". Error: " << (direct - half_left - half_right) / direct << std::endl;
            return direct;
        }

        double half_left = integrate(f, a, (a+b)/2.0, U());
        double half_right = integrate(f, (a+b)/2.0, b, U());

        if ( fabs(direct) > 1e-8 )
        {
            //check relative tolerance:
            if ( fabs( (direct - half_left - half_right) / direct) > 1e-8 )
            {
                //relative tolerance too large - split integration interval:
                return integrate(f, a, (a+b)/2.0, adaptiverule, recursion_depth + 1) + integrate(f, (a+b)/2.0, b, adaptiverule, recursion_depth + 1);
            }

            log::warn<log_integrate>()  << "Precision reached: a: " << a << ", b: " << b << std::endl;
            return half_left + half_right;
        }

        if ( fabs(direct - (half_left + half_right)) > 1e-8 )
        {
            //absolute tolerance too large - split integration interval:
            return integrate(f, a, (a+b)/2.0, adaptiverule, recursion_depth + 1) + integrate(f, (a+b)/2.0, b, adaptiverule, recursion_depth + 1);
        }

        return half_left + half_right;
    }

    /** @brief A binder for reducing an n-ary function to an n-1-ary function by integration over one variable */
    template <typename T, typename IntRule>
    class IntegrationBinder
    {
        public:
            IntegrationBinder(T const & func, double a, double b) : func_(func), a_(a), b_(b) {};

            double operator()(double x) const
            {
                BindFirst<T> integrand(func_, x);
                return integrate(integrand, a_, b_, IntRule());
            }


        private:
            T const & func_;
            double a_;
            double b_;
    };



    /** @brief Convenience overload for two-dimensional integration
     *
     * \\int_{a1}^{b1} \\int_{a2}^{b2} func2integrate(x,y) dy dx
     *
     * @tparam IntRuleTag     The integration rule to be used
     * @param  a1             Lower integration bound for the first variable
     * @param  b1             Upper integration bound for the second variable
     * @param  a2             Lower integration bound for the first variable
     * @param  b2             Upper integration bound for the second variable
     * @param  func2integrate The integrand
     */
    template <typename T, typename IntRuleTag>
    double integrate2D(double a1, double b1,
                       double a2, double b2, T const & func2integrate, IntRuleTag)
    {
        IntegrationBinder<T, IntRuleTag> innerIntegral(func2integrate, a2, b2);

        double result = integrate(innerIntegral, a1, b1, IntRuleTag());
        log::debug<log_integrate>()  << "check 0.25: " << innerIntegral(0.25) << std::endl;
        log::debug<log_integrate>()  << "check 0.5: " << innerIntegral(0.5) << std::endl;
        log::debug<log_integrate>()  << "check 0.75: " << innerIntegral(0.75) << std::endl;

        return result;
    }

  } //namespace math
} //namespace viennashe
#endif
