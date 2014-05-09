#ifndef VIENNASHE_VERSION_HPP
#define VIENNASHE_VERSION_HPP

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

#include <string>
#include <sstream>

/** @file  viennashe/version.hpp
    @brief Convenience functions for returning the current version of ViennaSHE
*/

namespace viennashe
{

  /** @brief Returns the current ViennaSHE major version as std::size_t   */
  inline std::size_t major_version()
  {
    return 1;
  }

  /** @brief Returns the current ViennaSHE minor version as std::size_t   */
  inline std::size_t minor_version()
  {
    return 2;
  }

  /** @brief Returns the current ViennaSHE revision number as std::size_t   */
  inline std::size_t revision_number()
  {
    return 0;
  }

  /** @brief Specifies whether the simulator is calibrated for some particular purpose.
   *
   * Note that a simulator requires calibration in order to yield predictive results!
   */
  inline bool calibrated() { return false; }

  /** @brief Returns the current ViennaSHE version string   */
  inline std::string version()
  {
    std::stringstream ss;
    ss << viennashe::major_version() << "."
       << viennashe::minor_version() << "."
       << viennashe::revision_number() << " " << ( (viennashe::calibrated() == true) ? "(calibrated)  " : "(uncalibrated)");
    return ss.str();
  }

  /** @brief  Prints the ViennaSHE preamble (header). Used in all the examples as well as the standalone-application */
  inline std::string preamble()
  {
    std::stringstream ss;
    ss << "*********************************************************" << std::endl;
    ss << "*             ViennaSHE " << viennashe::version() << "            *" << std::endl;
    ss << "*   A free open-source deterministic Boltzmann solver   *" << std::endl;
    ss << "*                      provided by                      *" << std::endl;
    ss << "*             Institute for Microelectronics            *" << std::endl;
    ss << "*    Institute for Analysis and Scientific Computing    *" << std::endl;
    ss << "*             Vienna University of Technology           *" << std::endl;
    ss << "*            http://viennashe.sourceforge.net/          *" << std::endl;
    ss << "*********************************************************" << std::endl;
    return ss.str();
  }
}


#endif

