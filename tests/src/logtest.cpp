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

/** \file logtest.cpp Contains a compilation test for the logging environment
 *  \test A simple, "Does it compile?"-test for the viennashe::log namespace
 */



#include <cstdlib>
#include <iostream>

#include "viennashe/log/log.hpp"

using namespace viennashe;

struct my_key_enabled
{
  enum { enabled = true };
};

struct my_key_disabled
{
  enum { enabled = false };
};

int main()
{
  std::cout << " ** BEGIN OF TEST **" << std::endl << std::endl;

  std::cout << "LOGTEST ... this is per stdout using fprintf ..."  << std::endl;
  std::cout << "Testing error, warn, info, debug now: " << std::endl << std::endl;

  log::error() << "ERROR VISIBLE " << std::endl;
  log::warn()  << "WARN  VISIBLE " << std::endl;
  log::info()  << "INFO  VISIBLE " << std::endl;
  log::debug() << "DEBUG NOT VISIBLE " << std::endl;

  log::error<my_key_enabled>() << "KEY: ERROR VISIBLE " << std::endl;
  log::warn<my_key_enabled>()  << "KEY: WARN  VISIBLE " << std::endl;
  log::info<my_key_enabled>()  << "KEY: INFO  VISIBLE " << std::endl;
  log::debug<my_key_enabled>() << "KEY: DEBUG NOT VISIBLE " << std::endl;

  log::error<my_key_disabled>() << "KEY: ERROR NOT VISIBLE " << std::endl;
  log::warn<my_key_disabled>()  << "KEY: WARN  NOT VISIBLE " << std::endl;
  log::info<my_key_disabled>()  << "KEY: INFO  NOT VISIBLE " << std::endl;
  log::debug<my_key_disabled>() << "KEY: DEBUG NOT VISIBLE " << std::endl;

  std::cout << std::endl;
  std::cout << " ** END OF TEST **" << std::endl;

  return EXIT_SUCCESS;
}

