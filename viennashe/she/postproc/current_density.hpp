#ifndef VIENNASHE_SHE_CURRENT_DENSITY_HPP
#define VIENNASHE_SHE_CURRENT_DENSITY_HPP

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



#include <math.h>
#include <vector>

#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/math/integrator.hpp"
#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/util/dual_box_flux.hpp"
#include "viennashe/postproc/current_density.hpp"
#include "viennashe/she/postproc/macroscopic.hpp"
#include "viennashe/she/assemble_common.hpp"

/** @file viennashe/she/postproc/current_density.hpp
    @brief Provides an accessor for the current density
*/

namespace viennashe
{
  namespace she
  {

    namespace detail
    {

      template <typename DeviceT, typename SHEQuantityT>
      class current_on_facet_by_ref_calculator
      {
          typedef viennashe::math::sparse_matrix<double>   CouplingMatrixType;

        public:
          current_on_facet_by_ref_calculator(DeviceT const & d,
                                             viennashe::config const & conf,
                                             SHEQuantityT const & quan)
           : device_(d), conf_(conf), quan_(quan),
             a_x(4, 4), a_y(4, 4), a_z(4, 4),
             b_x(4, 4), b_y(4, 4), b_z(4, 4)
          {
            fill_coupling_matrices(a_x, a_y, a_z,
                                   b_x, b_y, b_z,
                                   1);
          }

          std::vector<double> operator()(viennagrid_element_id facet) const
          {
            viennagrid_dimension cell_dim = viennagrid_topological_dimension_from_element_id(facet) + 1;

            viennashe::math::SphericalHarmonic Y_00(0,0);
            typename viennashe::config::dispersion_relation_type dispersion = conf_.dispersion_relation(quan_.get_carrier_type_id());

            viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device_.mesh(), facet, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

            std::vector<double> ret(3);

            if (cells_on_facet_begin + 1 == cells_on_facet_end)
              return ret;

            viennagrid_element_id c1 = cells_on_facet_begin[0];
            viennagrid_element_id c2 = cells_on_facet_begin[1];

            std::vector<viennagrid_numeric> normal_vector(3);
            viennagrid_numeric centroid_1[3] = { 0 };
            viennagrid_numeric centroid_2[3] = { 0 };
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c1, centroid_1));
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(device_.mesh(), c2, centroid_2));
            normal_vector[0] = centroid_2[0] - centroid_1[0];
            normal_vector[1] = centroid_2[1] - centroid_1[1];
            normal_vector[2] = centroid_2[2] - centroid_1[2];

            viennagrid_numeric norm;
            viennagrid_norm_2(3, &(normal_vector[0]), &norm);
            normal_vector[0] /= norm;
            normal_vector[1] /= norm;
            normal_vector[2] /= norm;

            double velocity_in_normal = 0;
            double polarity = (quan_.get_carrier_type_id() == ELECTRON_TYPE_ID) ? -1.0 : 1.0;

            double m_t = viennashe::materials::si::transverse_effective_mass(quan_.get_carrier_type_id());
            double m_d = viennashe::materials::si::dos_effective_mass(quan_.get_carrier_type_id());
            double m_l = viennashe::materials::si::longitudinal_effective_mass(quan_.get_carrier_type_id());

            for (std::size_t index_H = 1; index_H < quan_.get_value_H_size()-1; ++index_H)
            {
              long index = quan_.get_unknown_index(facet, index_H);
              if (index < 0)
                continue;

              double factor = 0;
              switch (conf_.she_discretization_type())
              {
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
                factor = integral_vZ(quan_, dispersion, facet, index_H) * dispersion.symmetry_factor();
                break;
              case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
                factor = integral_v(quan_, dispersion, facet, index_H) * dispersion.symmetry_factor();
                break;
              default: throw std::runtime_error("carrier_density_wrapper_by_reference::operator(): Unknown SHE discretization type!");
              }

              velocity_in_normal += (  a_x(0, 1) * quan_.get_values(facet, index_H)[0]
                                     + a_x(0, 2) * quan_.get_values(facet, index_H)[1]
                                     + a_x(0, 3) * quan_.get_values(facet, index_H)[2] ) * factor * sqrt(m_d / m_t) * normal_vector[0];

              if (cell_dim > 1)
                velocity_in_normal += (  a_y(0, 1) * quan_.get_values(facet, index_H)[0]
                                       + a_y(0, 2) * quan_.get_values(facet, index_H)[1]
                                       + a_y(0, 3) * quan_.get_values(facet, index_H)[2] ) * factor * sqrt(m_d / m_t) * normal_vector[1];

              if (cell_dim > 2)
                velocity_in_normal += (  a_z(0, 1) * quan_.get_values(facet, index_H)[0]
                                       + a_z(0, 2) * quan_.get_values(facet, index_H)[1]
                                       + a_z(0, 3) * quan_.get_values(facet, index_H)[2] ) * factor * sqrt(m_d / m_l) * normal_vector[2];
            }

            ret[0] = -polarity * viennashe::physics::constants::q * velocity_in_normal / Y_00(0,0);
            return ret;
          }


          SHEQuantityT const & get_quan() const { return this->quan_; }


        private:
          DeviceT      const & device_;
          viennashe::config conf_;
          SHEQuantityT const & quan_;

          CouplingMatrixType a_x;
          CouplingMatrixType a_y;
          CouplingMatrixType a_z;

          CouplingMatrixType b_x;
          CouplingMatrixType b_y;
          CouplingMatrixType b_z;
      };

      /** @brief A simple functor to hold a single value. Operator() should not take *any* arguments */
      template <typename ValueType>
      class value_holder_functor
      {
        public:
          typedef ValueType   value_type;

          template <typename T>
          void operator()(T const &, value_type val) const { value_ = val; }

          value_type operator()() const { return value_; }

        private:
          mutable value_type value_;
      };

    } // namespace detail


    /** @brief Accessor class providing the carrier velocity inside the device */
    template <typename DeviceType,
              typename SHEQuantity>
    class current_density_wrapper
    {
      public:
        typedef std::vector<double>       value_type;

        current_density_wrapper(DeviceType const & device,
                                viennashe::config const & conf,
                                SHEQuantity const & quan)
          : device_(device), quan_(quan), facet_evaluator_(device, conf, quan)
        { }

        current_density_wrapper(current_density_wrapper const & o)
          : device_(o.device_), quan_(o.quan_), facet_evaluator_(o.facet_evaluator_) {}

        /** @brief Functor interface returning the current density magnitude on the facet */
        std::vector<double> operator()(viennagrid_element_id cell_or_facet) const
        {
          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device_.mesh(), &cell_dim));

          std::vector<double> ret(3);
          if (viennagrid_topological_dimension_from_element_id(cell_or_facet) == cell_dim) // cell
          {
            if (!viennashe::materials::is_semiconductor(device_.get_material(cell_or_facet)))
              return ret;

            viennashe::util::value_holder_functor<std::vector<double> > result;

            viennashe::util::dual_box_flux_to_cell(device_,
                                                   cell_or_facet,
                                                   result, facet_evaluator_);
            return result();

          }
          else if (viennagrid_topological_dimension_from_element_id(cell_or_facet) == cell_dim - 1) // facet
          {
            ret[0] = facet_evaluator_(cell_or_facet)[0];
          }
          else
            throw std::runtime_error("current_density_wrapper::operator(): invalid element dimension!");

          return ret;
        }

      private:
        DeviceType const & device_;
        SHEQuantity  quan_;
        detail::current_on_facet_by_ref_calculator<DeviceType, SHEQuantity> facet_evaluator_; // REFERENCES quan_ and dispersion_
    };


    /** @brief Convenience function for writing the current density to the container provided.
     *
     * @param device           The device (includes a ViennaGrid mesh) on which simulation is carried out
     * @param quan             The SHE quantities. Passed for compatibility with other write_XYZ_to_domain() routines.
     * @param conf             The simulator configuration
     * @param container        Container for the current density vector
     */
    template <typename DeviceType, typename SHEQuantity>
    void write_current_density_to_quantity_field(DeviceType const & device,
                                                 viennashe::config const & conf,
                                                 SHEQuantity const & quan,
                                                 viennagrid_quantity_field field)
    {
      current_density_wrapper<DeviceType, SHEQuantity> current_wrapper(device, conf, quan);

      viennashe::write_macroscopic_quantity_to_quantity_field(device, current_wrapper, field);
    }


    /**
     * @brief Checks current conservation for SHE. Writes information using log::info().
     * @param device The device
     * @param conf   The simulator configuration
     * @param quan   The SHE quantities
     */
    template <typename DeviceT, typename SHEQuantity>
    void check_current_conservation(DeviceT const & device,
                                    viennashe::config const & conf,
                                    SHEQuantity const & quan)
    {
      detail::current_on_facet_by_ref_calculator<DeviceT, SHEQuantity> current_on_facet(device, conf, quan);

      viennashe::check_current_conservation(device, current_on_facet);
    }
  } //namespace she
} //namespace viennashe

#endif
