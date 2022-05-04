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
#include "libviennashe/include/material.h"


#ifdef	__cplusplus
extern "C" {
#endif

/** @brief  Returns the ID of silicon */
viennasheErrorCode viennashe_get_silicon_id (viennashe_material_id * id)
{
  if(id != NULL) *id = viennashe::materials::si::id;
  else return 1;
  return 0;
}

/** @brief Returns the ID of ideal conductors */
viennasheErrorCode viennashe_get_metal_id   (viennashe_material_id * id)
{
  if(id != NULL) *id = viennashe::materials::metal::id;
  else return 1;
  return 0;
}

/** @brief Returns the ID of SiO2 */
viennasheErrorCode viennashe_get_sio2_id    (viennashe_material_id * id)
{
  if(id != NULL) *id = viennashe::materials::sio2::id;
  else return 1;
  return 0;
}

/** @brief Returns the ID of HfO2 */
viennasheErrorCode viennashe_get_hfo2_id    (viennashe_material_id * id)
{
  if(id != NULL) *id = viennashe::materials::hfo2::id;
  else return 1;
  return 0;
}

#ifdef	__cplusplus
}
#endif

