#ifndef VIENNASHE_SHE_RESCALING_HPP
#define VIENNASHE_SHE_RESCALING_HPP

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
#include <iostream>
#include <math.h>
#include <cfloat>

// viennashe
#include "viennashe/config.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/she/she_quantity.hpp"

// viennagrid
#include "viennagrid/mesh/mesh.hpp"

/** @file viennashe/she/rescaling.hpp
    @brief Rescales the linear system such that unknowns are approximately the same order of magnitude.
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Initializes the unknown scaling. Returns a vector holding the scaling factors
     *
     * @param device           The device (includes a ViennaGrid mesh) on which a simulation is carried out
     * @param conf             The simulator config
     * @param quan             The SHE quantities
     * @param scaling_vector   The vector to be filled with the scaling factors
     */
    template <typename DeviceType,
              typename SHEQuantity,
              typename VectorType>
    void setup_unknown_scaling(DeviceType const & device,
                               viennashe::config const & conf,
                               SHEQuantity & quan,
                               VectorType & scaling_vector)
    {
      typedef typename DeviceType::mesh_type                                          MeshType;
      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type        CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type           CellIterator;

      CellContainer cells(device.mesh());
      for (CellIterator cit = cells.begin();
          cit != cells.end();
          ++cit)
      {
        double T = device.get_lattice_temperature(*cit);

        for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
        {
          long index_base = quan.get_unknown_index(*cit, index_H);
          if (index_base < 0)
            continue;

          long num_per_node = static_cast<long>(quan.get_unknown_num(*cit, index_H));
          for (long i=0; i<num_per_node; ++i)
          {
            if (conf.she_scaling_type() == SHE_SCALING_KINETIC_ENERGY)
              scaling_vector[std::size_t(index_base + i)] = exp(- quan.get_kinetic_energy(*cit, index_H) / (viennashe::physics::constants::kB * T));
            else
              scaling_vector[std::size_t(index_base + i)] = exp(- quan.get_value_H(index_H) / (viennashe::physics::constants::kB * T));
          }
        }
      }
    }

    /** @brief Scales the system matrix (right-preconditioner) with respect to the rescaled unknowns
     *
     * @param system_matrix       The system matrix to be scaled with a diagonal right-preconditioner
     * @param scaling             The values of the diagonal of the right-preconditioner
     */
    template <typename MatrixType,
              typename VectorType>
    void rescale_system(MatrixType & system_matrix,
                        VectorType const & scaling)
    {
      typedef typename MatrixType::row_type      RowType;
      typedef typename MatrixType::iterator2     AlongRowIterator;

      for (std::size_t i=0; i<system_matrix.size1(); ++i)
      {
        RowType & row_i = system_matrix.row(i);
        for (AlongRowIterator iter  = row_i.begin();
                              iter != row_i.end();
                            ++iter)
        {
          iter->second *= scaling[iter->first];
        }
      }
    }


  } //namespace she
} //namespace viennashe

#endif
