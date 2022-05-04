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

// C++ includes
#include "libviennashe/src/viennashe_all.hpp"

// C includes
#include "libviennashe/include/error.h"
#include "libviennashe/include/libviennashe.h"


#ifdef	__cplusplus
extern "C" {
#endif

/**
 *
 * @return
 */
viennasheErrorCode viennashe_initalize(void)
{
  try
  {
    viennashe::log::info() << viennashe::preamble() << std::endl;
    //
    // If you are using MPI: Initialize MPI here!
    //

    return 0;
  }
  catch(...)
  {
    return -1;
  }
}

/**
 *
 * @return
 */
viennasheErrorCode viennashe_finalize(void)
{
  //
  // If you are using MPI: Finalize MPI here!
  //

  viennashe::log::info() << std::endl;
  viennashe::log::info() << "*********************************************************" << std::endl;
  viennashe::log::info() << "*                ViennaSHE finished                     *" << std::endl;
  viennashe::log::info() << "*********************************************************" << std::endl;
  return 0;
}

/**
 * @brief If the ecode != 0, the routine prints error information onto the screen.
 * @param ecode The error code
 * @return The argument ecode
 */
viennasheErrorCode viennashe_error(viennasheErrorCode ecode)
{
  if (ecode != 0)
  {
    if(ecode > 0) viennashe::log::error() << "Found error in the last routine. Error is likely to be due to argument number " << ecode << std::endl;
    else          viennashe::log::error() << "Encountered a severe internal error in the last routine. " << std::endl;

    return ecode;
  }

  return 0;
}

#ifdef	__cplusplus
}
#endif
