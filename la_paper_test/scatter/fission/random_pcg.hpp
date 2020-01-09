/*
 *  Libre Monte Carlo Transport (LMCT)
 *
 *  Copyright (c) 2020 Hunter Belanger (hunter.belanger@gmail.com)
 *  All rights reserved.
 *
 *  Released under the BSD 3-Clause license. Please refer to the
 *  LICENSE file in the root directory of this distribution.
 *
 *  File: random_pcg.hpp
 *
 *  Contents: Random number generatior
 *
 * */
#ifndef LMCT_RANDOM_PCG_H
#define LMCT_RANDOM_PCG_H

#include"pcg_random.hpp"

#include<random>
#include<vector>

namespace lmct {
  
  // ---------------------------------------------------------------------------
  // PCG
  //   Wrapper class for the pcg32 random number generator. Stores the current
  //   seed of the engine, and also has class method to produce a random double
  //   over the interval [0,1).
  // ---------------------------------------------------------------------------
  class PCG {
    public:
      PCG() {}
      PCG(uint64_t _seed) {engine.set_seed(_seed);}
      ~PCG() = default;

      // -----------------------------------------------------------------------
      // seed(uint64_t _seed)
      //   Directly sets seed (state_) for the pcg32 engine.
      // -----------------------------------------------------------------------
      void seed(uint64_t _seed) {
        engine.set_seed(_seed);
      }
      
      // -----------------------------------------------------------------------
      // get_seed()
      //   Returns the current seed of the pcg32 engine;
      // -----------------------------------------------------------------------
      uint64_t get_seed() {
        return engine.get_seed();
      }
      
      // -----------------------------------------------------------------------
      // rand()
      //   This advances the pcg32 generator and produces a double
      //   over the interval [0,1).
      // -----------------------------------------------------------------------
      double rand() {
        return unit_dist(engine);
      }

      // -----------------------------------------------------------------------
      // uniform(double a, double b)
      //   Returns a random double with a uniform distribution over the
      //   interval [a, b)
      //
      //   f(x|a,b) = 1 / (b - a)
      // -----------------------------------------------------------------------
      double uniform(double a, double b) {
        std::uniform_real_distribution<double> dist(a,b);
        return dist(engine);
      }

      // -----------------------------------------------------------------------
      // normal(double mu, double sigma)
      //   Returns a random double from the normal (Gaussian) distribution,
      //   defined by average value mu, and std sigma
      //
      //   f(x|mu,sigma) = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x - mu)/sigma)^2)
      // -----------------------------------------------------------------------
      double normal(double mu, double sigma) {
        std::normal_distribution<double> dist(mu, sigma);
        return dist(engine);
      }

      // -----------------------------------------------------------------------
      // exponential(double lambda)
      //   Returns a random double from the exponential distribution,
      //   defined by the average value 1/lambda.
      //
      //   f(x|lambda) = lambda * exp(-lambda * x)
      // -----------------------------------------------------------------------
      double exponential(double lambda) {
        std::exponential_distribution<double> dist(lambda);
        return dist(engine);
      }

      // -----------------------------------------------------------------------
      // discrete(std::vector<double> weights)
      //   Returns an integer over range [0,weights.size() - 1]
      //   where probability of each integer is defined by the weight
      //
      //   P(i|w_0, w_1, ... w_k) = w_i / Sum[j = 1 to k](w_j)
      // -----------------------------------------------------------------------
      int discrete(std::vector<double> weights) {
        std::discrete_distribution<int> dist(weights.begin(), weights.end());
        return dist(engine);
      }

    private:
      pcg32 engine;

      std::uniform_real_distribution<double> unit_dist;

  }; // PCG

} // lmct

#endif // LMCT_RANDOM_PCG_H
