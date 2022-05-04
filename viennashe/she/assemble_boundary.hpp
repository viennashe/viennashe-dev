#ifndef VIENNASHE_SHE_ASSEMBLE_BOUNDARY_HPP
#define VIENNASHE_SHE_ASSEMBLE_BOUNDARY_HPP

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


// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/norm.hpp"
#include "viennagrid/algorithm/volume.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// viennashe
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/config.hpp"
#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/assemble_common.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/util/filter.hpp"

#include "viennashe/log/log.hpp"

/** @file viennashe/she/assemble_streaming.hpp
    @brief Assembly of the free-streaming operator is implemented here
*/

namespace viennashe
{
  namespace she
  {

    /** @brief Worker function for the assembly of the free streaming operator. Handles the assembly of both even and odd unknowns.
     *
     * @tparam DeviceType        The device descriptor class
     * @tparam SHEQuantity       The SHE quantity type giving access to respective data (unknown indices, etc.)
     * @tparam MatrixType        The system matrix type
     * @tparam VectorType        The load vector type
     * @tparam CellType          The topological cell element type for boundary assembly
     * @tparam CouplingMatrixType Type of the coupling matrices a_{l,m}^{l',m'} and b_{l,m}^{l',m'}
     */
    template <typename DeviceType,
              typename SHEQuantity,
              typename MatrixType,
              typename VectorType,
              typename CellType,
              typename CouplingMatrixType>
    void assemble_boundary_on_box(DeviceType const & device,
                                  viennashe::config const & conf,
                                  SHEQuantity const & quan,
                                  MatrixType & A,
                                  VectorType & b,
                                  CellType const & cell, std::size_t index_H,
                                  CouplingMatrixType const & coupling_identity)
    {
      const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

      long row_index = quan.get_unknown_index(cell, index_H);

      if (row_index < 0) //no unknown here
        return;

      long expansion_order = static_cast<long>(quan.get_expansion_order(cell,  index_H));

      if (device.has_contact_potential(cell)) // any boundary conditions required here?
      {
        typename viennashe::config::dispersion_relation_type dispersion = conf.dispersion_relation(quan.get_carrier_type_id());

        if (conf.she_boundary_conditions().type() == viennashe::BOUNDARY_DIRICHLET) // Dirichlet boundary conditions
        {
          double bnd_value = quan.get_boundary_value(cell, index_H);
          double height = box_height(quan, cell, index_H);

          double coefficient = 0;
          switch (conf.she_discretization_type())
          {
          case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
            coefficient = height * averaged_density_of_states(quan, dispersion, cell, index_H);
            break;
          case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
            coefficient = height;
            break;
          default: throw std::runtime_error("assemble_boundary_on_box(): Unknown SHE discretization type!");
          }

          viennashe::util::add_block_matrix(A, std::size_t(row_index), std::size_t(row_index),
                                            coefficient,
                                            coupling_identity,
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID)
                                        );

          write_boundary(b, std::size_t(row_index),
                          -bnd_value*coefficient,
                          coupling_identity,
                          viennashe::math::spherical_harmonics_iterator(1, viennashe::math::EVEN_HARMONICS_ITERATION_ID)  //equilibrium distribution is zeroth-order
                        );
          return;
        }
        else // generation/recombination boundary conditions
        {
          //
          // boundary term (for even unknowns):  dg/dt + L - Q + (g - g^eq) / tau = 0  as a volume contribution on RHS
          //

          double tau    = conf.she_boundary_conditions().generation_recombination_rate();
          double height = box_height(quan, cell, index_H);
          double bnd_value  = quan.get_boundary_value(cell, index_H);
          double box_volume = viennagrid::volume(cell);

          double coefficient = 0;
          switch (conf.she_discretization_type())
          {
          case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
            coefficient = box_volume * height / tau * averaged_density_of_states(quan, dispersion, cell, index_H);
            break;
          case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
            coefficient = box_volume * height / tau;
            break;
          default: throw std::runtime_error("assemble_boundary_on_box(): Unknown SHE discretization type!");
          }

          viennashe::util::add_block_matrix(A, std::size_t(row_index), std::size_t(row_index),
                                            coefficient,
                                            coupling_identity,
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                            viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID)
                                        );

          write_boundary(b, std::size_t(row_index),
                          - bnd_value * coefficient,  //Note: sign choice is the same as for add_block_matrix(), i.e. write_boundary() inverts the sign for the rhs contribution internally
                          coupling_identity,
                          viennashe::math::spherical_harmonics_iterator(1, viennashe::math::EVEN_HARMONICS_ITERATION_ID)  //equilibrium distribution is zeroth-order
                        );

          if (log_assemble_free_streaming_operator::debug)
          {
            log::debug<log_assemble_free_streaming_operator>() << "* assemble_boundary_on_box(): boundary contribution:" << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "      bnd_value * box_volume * height / tau" << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "  at (" << row_index << ")" << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "    bnd_value = " << bnd_value << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "   box_volume = " << box_volume << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "       height = " << height << std::endl;
            log::debug<log_assemble_free_streaming_operator>() << "          tau = " << tau << std::endl;
          }

          if (with_full_newton)
          {
            // Only contributions to the residual, because (f - f^eq)/tau does neither have an explicit dependence on the potential, nor does it spatially couple.
            viennashe::util::subtract_folded_block_vector(b,
                                                          std::size_t(row_index),
                                                          box_volume * height / tau,
                                                          coupling_identity,
                                                          viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                                          viennashe::math::spherical_harmonics_iterator(expansion_order, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                                          quan.get_values(cell, index_H));
          }
        }

      } // if boundary_terms

    } //assemble_boundary_on_box



  } //namespace she
} //namespace viennashe

#endif
