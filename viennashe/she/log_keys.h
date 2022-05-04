#ifndef VIENNASHE_SHE_LOG_KEYS_H
#define VIENNASHE_SHE_LOG_KEYS_H

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

/** @file viennashe/she/log_keys.h
    @brief Defines the log keys used within the viennashe::she namespace
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Configuration class for logging during the linear solver phase in a SHE simulation */
    struct log_she_solver { enum { enabled = true, debug = false }; };

    struct log_newton_she { enum { enabled = true }; };

    //
    // scattering:
    //
    /** @brief Configuration class for logging when assembling acoustic phonon scattering */
    struct log_acoustic_phonon_scattering { enum { enabled = false }; };

    /** @brief Configuration class for logging when assembling optical phonons scattering */
    struct log_optical_phonon_scattering { enum { enabled = false }; };

    /** @brief Configuration class for logging when assembling impact ionization scattering*/
    struct log_impact_ionization_scattering { enum { enabled = false }; };

    /** @brief Configuration class for logging when assembling electron-electron scattering*/
    struct log_assemble_ee_scattering { enum { enabled = true }; };

    /** @brief Configuration class for logging when assembling any scattering mechanisms */
    struct log_assemble_scattering_operator { enum { enabled = false }; };

    //
    // misc:
    //

    /** @brief Configuration class for logging when smoothing expansion orders for mapping */
    struct log_smooth_expansion_order { enum { enabled = true }; };

    /** @brief Configuration class for logging when distributing unknown indices over device (mapping) */
    struct log_mapping { enum { enabled = true }; };

    /** @brief Configuration class for logging when running the linear solver */
    struct log_linear_solver { enum { enabled = true, debug = false }; };

    /** @brief Configuration class for logging when filling coupling matrices */
    struct log_fill_coupling_matrices { enum { enabled = true }; };

    /** @brief Configuration class for logging when recovering odd unknowns from the computed even f_{l,m} */
    struct log_recover_odd_unknowns { enum { enabled = true }; };

    /** @brief Configuration class for logging when computing the coupling matrix in a particular direction */
    struct log_coupling_matrix_in_direction { enum { enabled = true }; };

    /** @brief Configuration class for logging when assembling the free streaming operator */
    struct log_assemble_free_streaming_operator { enum { enabled = true, debug = false }; };

    /** @brief Configuration class for logging when assembling the free streaming operator */
    struct log_assemble_all { enum { enabled = false, debug = false }; };

    /** @brief Configuration class for logging when using adaptive SHE */
    struct log_adaptive_she { enum { enabled = true }; };

    /** @brief Configuration class for logging when dealing with traps */
    struct log_traps { enum { enabled = false }; };

    /** @brief Configuration class for logging when dealing with traps */
    struct log_transfer_to_new_h_space { enum { enabled = false }; };

  }
} // namespace viennashe


#endif

