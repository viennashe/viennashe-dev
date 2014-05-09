#ifndef VIENNASHE_MATH_RANDOM_HPP
#define VIENNASHE_MATH_RANDOM_HPP

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
#include <cmath>
#include <stdio.h>
#include <stdint.h>

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/log/log.hpp"


/** @file viennashe/math/random.hpp
    @brief Implementation of pseudo-random number generators
*/

namespace viennashe
{
  namespace math
  {

    /** @brief The basic interface for pseudo-random number generators */
    template< typename ValueT = double>
    class rand_generator_base
    {
    public:
      virtual ~rand_generator_base() { }

      /**
       * @brief The basic random number generator interface (virtual abstract)
       * @return A pseudo-randum number between 0 and 1.0 (either 0 or 1 can be excluded!)
       */
      virtual ValueT operator()()  = 0;
    };

    /** @brief A pseudo-random number generator implementation based on std::rand() */
    template< typename ValueT = double>
    class std_rand_generator : public rand_generator_base<ValueT>
    {
    public:
        /**
         * @brief Sets the seed of the pseudo-random number generator (std::srand)
         * @param s The seed to get reproduceable sequences of numbers
         */
        void seed(unsigned int s)  { std::srand(s);  }

        /**
         * @brief Calls std::rand divides by RAND_MAX and casts to ValueT
         * @return A pseudo-random number in [0,1]
         */
        ValueT operator()()   { return static_cast<ValueT>(std::rand()*(1.0/RAND_MAX)); }
    };

    /** @brief A pseudo-random number generator implementation using the merseinne twister MT 19937
     *         cf. Makoto Matsumoto, Takuji Nishimura: Mersenne twister.
     *         "A 623-dimensionally equidistributed uniform pseudorandom number generator."
     *         ACM Transactions on Modeling and Computer Simulation
     *         8, 1998
     */
    template< typename ValueT = double>
    class merseinne_twister_generator
    {
    public:
      /**
       * @brief Default CTOR to ensure default constructability
       */
      merseinne_twister_generator() : index_(N_ + 1), seed_(5489UL) { this->init(); }

      /** @brief CTOR, which accepts the seed to the generator */
      merseinne_twister_generator(unsigned int seed) : index_(N_ + 1), seed_(seed) { this->init(); }

      /**
       * @brief Sets the seed of the pseudo-random number generator
       * @param s The seed to get reproduceable sequences of numbers
       */
      void seed(unsigned int s) { seed_ = s; }

      /** @brief Returns a pseudo-random number using MT19937.
       *         Not const, since the internal state is being changed!
       */
      ValueT operator()()
      {
        const uint32_t HI_ = 0x80000000;
        const uint32_t LO_ = 0x7fffffff;
        const uint32_t A_[2] = { 0, 0x9908b0df };

        uint32_t value = 0;

        // Init if not done yet
        if (index_ == N_+1) this->init();

        // Next state
        if (index_ >= N_)
        {
          uint32_t h = 0;
          uint32_t i = 0;

          for (i = 0; i < N_-M_; ++i)
          {
            h = (y_[i] & HI_) | (y_[i+1] & LO_);
            y_[i] = y_[i+M_] ^ (h >> 1) ^ A_[h & 1];
          }
          for ( ; i < N_-1; ++i)
          {
            h = (y_[i] & HI_) | (y_[i+1] & LO_);
            y_[i] = y_[i+(M_-N_)] ^ (h >> 1) ^ A_[h & 1];
          }

          h = (y_[N_-1] & HI_) | (y_[0] & LO_);
          y_[N_-1] = y_[M_-1] ^ (h >> 1) ^ A_[h & 1];
          index_ = 0;
        }
        value = y_[index_++];

        // Tempering
        value ^= (value >> 11);
        value ^= (value <<  7) & 0x9d2c5680;
        value ^= (value << 15) & 0xefc60000;
        value ^= (value >> 18);

        return (static_cast<ValueT>(value)) * (1.0/4294967296.0) /* * 1/2^32 */;
      }

    private:

      /** @brief Initializes the merseinne twister */
      void init()
      {
        y_[0] = seed_; // set seed
        for (uint32_t i = 1; i < N_; ++i)
          y_[i] = (1812433253UL * (y_[i-1] ^ (y_[i-1] >> 30)) + i);
      }

      static const uint32_t N_  = 624;
      static const uint32_t M_  = 397;

      uint32_t index_;
      uint32_t seed_;
      uint32_t y_[N_];
    };


    /** @brief A uniform distribution of random numbers */
    template <typename ResultT = int>
    class uniform_distribution
    {
    public:

      /**
       * @brief Returns a uniformly distributed pseudo-random number
       * @param gen A random number generator implementing the concept of rand_generator_base
       * @param from Minimum of the uniform distribution
       * @param to Maximum of the uniform distribution
       * @return Uniformly distributed random numbers
       */
      template<typename RandGeneratorType >
      ResultT operator()(RandGeneratorType & gen, ResultT from, ResultT to)
      {
          return static_cast<ResultT>(gen() * (to - from) + from);
      }

    }; // uniform_distribution


    /** @brief A poisson (exponential) distribution of random numbers */
    template <typename RealT = double, typename ResultT = int>
    class poisson_distribution
    {
    private:
      RealT mean_; // The mean of the distribution
      RealT exp_mean_; // Equals std::exp(- mean_); Used for caching
    public:

      poisson_distribution(RealT const & mean_value)
        : mean_(mean_value)
      {
        this->exp_mean_ = std::exp(- this->mean_);
      }

      RealT mean() const  { return this->mean_; }

      /**
       * @brief Returns a poisson distributed random number
       * @param gen A random number generator implementing the concept of rand_generator_base
       * @return Poisson distributed random number
       */
      template <typename RandGeneratorType>
      ResultT operator()(RandGeneratorType & gen)
      {
        RealT product = RealT(0);
        for(ResultT m = 0; ; ++m)
        {
          product += std::log(gen());
          if(product <= -this->mean_)
            return m;
        }
      }

    }; // poisson_distribution


    /*
     * @brief Pseudo-random number generator for normally distributed random numbers
     *        Using the Box-Muller transform
     *        Algorithm from: Numerical Recipes, The Art of Scientific Computing, Third Edition, Page 364/365
     */
    template <typename RealT = double, typename ResultT = double>
    class normal_distribution
    {
    private:
      RealT mean_; // The mean value
      RealT sigma_; // The standard deviation
      RealT storedvalue_; // For internal caching

    public:
      normal_distribution(RealT const & mean,
                          RealT const & sigma)
        : mean_(mean), sigma_(sigma), storedvalue_(0.0) {  }

      RealT mean()  const  { return this->mean_;  }
      RealT sigma() const  { return this->sigma_; }

      /**
       * @brief Returns a Gaussian distributed random number; Algorithm: Box-Muller
       * @param gen A random number generator implementing the concept of rand_generator_base
       * @return A Gaussian distributed random number
       */
      template<typename RandGeneratorType>
      ResultT operator()(RandGeneratorType &gen)
      {
        double v1, v2, rsq, fac;
        if (storedvalue_ == 0.0)                    // We don't have an extra deviate handy, so ...
        {
          do
          {
            v1 = 2 * gen() - 1.0;                 // pick two uniform numbers in the square extending
            v2 = 2 * gen() - 1.0;                 //   from -1 to +1 in each direction
            rsq = v1*v1 + v2*v2;                  // see if they are in the unit circle,
          } while (rsq >= 1.0 || rsq == 0.0);      //    or try again.
          fac = std::sqrt(-2.0 * std::log(rsq) / rsq);// Now make the Box-Muller transformation to
          storedvalue_ = v1 * fac;                    //    get two normal deviates. Return one and
          return this->mean_ + this->sigma_*v2*fac;   //    safe the other for next time.
        }
        else                                        // We have an extra deviate handy ...
        {
          fac = storedvalue_;
          storedvalue_ = 0.0; // reset
          return this->mean_ + this->sigma_*fac;     // ... so return it.
        }
      }

    }; // normal_distribution



  } // namespace math

} // namespace viennashe

#endif /* VIENNASHE_MATH_RANDOM_HPP */

