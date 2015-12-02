#ifndef VIENNASHE_SHE_ASSEMBLE_COMMON_HPP
#define VIENNASHE_SHE_ASSEMBLE_COMMON_HPP

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


// viennagrid
#include "viennagrid/viennagrid.h"

// viennashe
#include "viennashe/config.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/harmonics_coupling.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/misc.hpp"

#include "viennashe/she/exception.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"


/** @file viennashe/she/assemble_common.hpp
    @brief Common routines for the assembly of SHE equations.
*/

namespace viennashe
{
  namespace she
  {
    namespace detail
    {
      template <typename PotentialAccessorType,
                typename ConfigType, typename TagType>
      std::vector<long> get_potential_indices(PotentialAccessorType const & potential_index,
                                              viennagrid_mesh mesh,
                                              viennagrid_element_id elem)
      {
        std::vector<long> ret;

        if (viennagrid_topological_dimension_from_element_id(elem) > 0)
        {

          viennagrid_element_id *vertices_begin, *vertices_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, 0, &vertices_begin, &vertices_end));

          for (viennagrid_element_id *it=vertices_begin; it != vertices_end; ++it)
          {
            if (potential_index(viennagrid_index_from_element_id(*it)) >= 0)
              ret.push_back(potential_index(*it));
          }
        }
        else
        {
          if (potential_index(viennagrid_index_from_element_id(elem)) >= 0)
            ret.push_back(potential_index(viennagrid_index_from_element_id(elem)));
        }
        return ret;
      }



      /** @brief Static dispatcher for computing the connection between two cells adjacent to a facet. */
      template <typename MeshT>
      double cell_connection_length(MeshT const & mesh, viennagrid_element_id facet)
      {
        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

        viennagrid_element_id *cells_begin, *cells_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh, facet, cell_dim, &cells_begin, &cells_end));

        if (cells_end == cells_begin + 1)
        {
          assert(bool("Logic error: cell_connection_length() called for facet on boundary!"));
        }

        std::vector<double> centroid1(3), centroid2(3);
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, cells_begin[0], &(centroid1[0])));
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, cells_begin[1], &(centroid2[0])));

        viennagrid_dimension geo_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_geometric_dimension_get(mesh, &geo_dim));

        // compute 2-norm of distance between centroids
        viennagrid_numeric distance;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(geo_dim, &(centroid1[0]), &(centroid2[0]), &distance));

        return distance;
      }

      template <typename DeviceType>
      bool has_contact_potential(DeviceType const & device, viennagrid_element_id c)
      {
        return device.has_contact_potential(c);
      }

    }


    /** @brief Computes \\int_a^b Z/(hbar * k) dE, where a and b are the energy boundaries, Z is the density of states and hbar * k is the momentum.
     *
     * Uses the formula Z/(hbar * k) = 0.5 d(Zv)/dE, thus Z and v are evaluated at the energy boundaries a and b only
     *
     * @param energy_minus   Lower integration bound (might be smaller than zero -> set to zero)
     * @param energy_plus    Upper integration bound (might be smaller than zero -> set to zero)
     * @param dispersion     The dispersion relation used for the SHE simulation
     */
    template <typename DispersionRelation>
    double integral_Z_over_hk(double energy_minus, double energy_plus, DispersionRelation const & dispersion)
    {
      if (energy_plus < energy_minus)
      {
        std::stringstream ss;
        ss << "Invalid interval in integral_Z_over_hk(): [" << energy_minus << ", " << energy_plus << "]";
        throw negative_integration_interval_length_exception(ss.str());
      }

      return 0.5 * ( dispersion.density_of_states(energy_plus) * dispersion.velocity(energy_plus)
                     - dispersion.density_of_states(energy_minus) * dispersion.velocity(energy_minus) );
    }


    /** @brief Returns dot(M, n), where M=(m1, m2, m3) is the vector of coupling matrices, and n = v2 - v1 is the directional vector (i.e. normal vector on box). One spatial dimension.
     *
     * @param m1   Coupling matrix in x-direction
     * @param m2   Coupling matrix in y-direction
     * @param m3   Coupling matrix in z-direction
     * @param c1   First  cell (no Delaunay criterion required!)
     * @param c2   Second cell (no Delaunay criterion required!)
     * @param polarity  A polarity tag for distinguishing between electrons and holes
     */
    template <typename MatrixType, typename PolarityTag>
    MatrixType coupling_matrix_in_direction_impl(MatrixType const & m1, MatrixType const & m2, MatrixType const & m3,
                                                 viennagrid_mesh mesh, viennagrid_element_id c1, viennagrid_element_id c2,
                                                 PolarityTag const & polarity)
    {
      std::vector<double> centroid1(3), centroid2(3);
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, c1, &(centroid1[0])));
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, c2, &(centroid2[0])));

      double norm = 0;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_distance_2(3, &(centroid1[0]), &(centroid2[0]), &norm));

      double alpha_1 = (centroid2[0] - centroid1[0]) / norm;
      double alpha_2 = (centroid2[1] - centroid1[1]) / norm;
      double alpha_3 = (centroid2[2] - centroid1[2]) / norm;

      //log::debug<log_coupling_matrix_in_direction>() << "Coefficients: " << p3;
      double m_t = viennashe::materials::si::transverse_effective_mass(polarity);
      double m_l = viennashe::materials::si::longitudinal_effective_mass(polarity);
      double m_d = viennashe::materials::si::dos_effective_mass(polarity);

      MatrixType result = m1 * (alpha_1 * sqrt(m_d / m_t));
      result += m2 * (alpha_2 * sqrt(m_d / m_t));
      result += m3 * (alpha_3 * sqrt(m_d / m_l));
      return result;
    }


    /** @brief Returns dot(M, n), where M=(m1, m2, m3) is the vector of coupling matrices, and n = v2 - v1 is the directional vector (i.e. normal vector on box)
     *
     * @param m1   Coupling matrix in x-direction
     * @param m2   Coupling matrix in y-direction
     * @param m3   Coupling matrix in z-direction
     * @param v1   First vertex (of Delaunay triangulation connecting two Voronoi boxes)
     * @param v2   Second vertex (of Delaunay triangulation connecting two Voronoi boxes)
     * @param polarity  A polarity tag for distinguishing between electrons and holes
     */
    template <typename MatrixType, typename PolarityTag>
    MatrixType coupling_matrix_in_direction(MatrixType const & m1, MatrixType const & m2, MatrixType const & m3,
                                            viennagrid_mesh mesh, viennagrid_element_id v1, viennagrid_element_id v2,
                                            PolarityTag const & polarity)
    {
      //typedef typename viennagrid::result_of::point<VertexType>::type             PointType;
      //typedef typename viennagrid::result_of::coordinate_system<PointType>::type  CoordinateSystemTag;
      //return coupling_matrix_in_direction_impl(m1, m2, m3, v1, v2, polarity, CoordinateSystemTag());  //using ViennaGrid vertices here
      return coupling_matrix_in_direction_impl(m1, m2, m3, mesh, v1, v2, polarity);  //using ViennaGrid vertices here
    }


    /** @brief Write Dirichlet boundary conditions to right hand side
     *
     * @param rhs              Load vector
     * @param row_start        Load vector index to start
     * @param prefactor        Scalar prefactor for 'coupling_matrix'
     * @param coupling_matrix  SHE coupling matrix
     * @param row_iter         Iterator over row indices in 'coupling_matrix'
     */
    template <typename MatrixType, typename VectorType, typename IteratorType>
    void write_boundary(VectorType & rhs,
                        std::size_t row_start,
                        double prefactor,
                        MatrixType const & coupling_matrix,
                        IteratorType row_iter)
    {
      assert(prefactor == prefactor && bool("Writing nan to boundary!"));

      //this function is called only for odd rows (more precisely: odd spherical harmonics) and even cols (more precisely: even spherical harmonics)
      //thus: lowest even-order harmonic is unequal to zero and goes to right hand side
      std::size_t row = row_start;
      for (; row_iter.valid(); ++row_iter)
      {
        rhs[row] -= prefactor * coupling_matrix(std::size_t(*row_iter), 0);
        ++row;
      }

    }



    /** @brief Returns the lower kinetic energy for the discretization box associated with the edge 'edge'. */
    template <typename DeviceType, typename SHEQuantity>
    double lower_kinetic_energy_at_facet(DeviceType const & device,
                                         SHEQuantity const & quan,
                                         viennagrid_element_id facet,
                                         std::size_t index_H)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), facet, cell_dim, &cells_begin, &cells_end));

      if (cells_begin + 1 == cells_end)
        return quan.get_kinetic_energy(cells_begin[0], index_H);

      if (index_H == 0)
        return quan.get_kinetic_energy(cells_begin[0], index_H) + quan.get_kinetic_energy(cells_begin[1], index_H) / 2.0;

      return (  quan.get_kinetic_energy(cells_begin[0], index_H - 1) + quan.get_kinetic_energy(cells_begin[1], index_H - 1)
              + quan.get_kinetic_energy(cells_begin[0], index_H)     + quan.get_kinetic_energy(cells_begin[1], index_H)
              ) / 4.0;
    }


    /** @brief Returns the upper kinetic energy of the discretization box for the edge 'edge' */
    template <typename DeviceType, typename SHEQuantity>
    double upper_kinetic_energy_at_facet(DeviceType const & device,
                                         SHEQuantity const & quan,
                                         viennagrid_element_id facet,
                                         std::size_t index_H)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), facet, cell_dim, &cells_begin, &cells_end));

      if (cells_begin + 1 == cells_end)
        return quan.get_kinetic_energy(cells_begin[0], index_H);

      if (index_H == quan.get_value_H_size() - 1)
        return quan.get_kinetic_energy(cells_begin[0], index_H) + quan.get_kinetic_energy(cells_begin[1], index_H) / 2.0;

      return (  quan.get_kinetic_energy(cells_begin[0], index_H + 1) + quan.get_kinetic_energy(cells_begin[1], index_H + 1)
              + quan.get_kinetic_energy(cells_begin[0], index_H)     + quan.get_kinetic_energy(cells_begin[1], index_H)
              ) / 4.0;
    }

    class force_on_facet_accessor
   {
    public:

      template <typename SHEQuantity>
      double operator()(SHEQuantity const & quan,
                        viennagrid_mesh mesh,
                        viennagrid_element_id const & c1,
                        viennagrid_element_id const & c2,
                        std::size_t index_H) const
      {
        std::vector<double> centroid1(3), centroid2(3);
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, c1, &(centroid1[0])));
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_centroid(mesh, c2, &(centroid2[0])));

        for (std::size_t i=0; i<centroid2.size(); ++i)
          centroid2[i] -= centroid1[i];

        double distance = 0;
        for (std::size_t i=0; i<centroid2.size(); ++i)
          distance += centroid2[i] * centroid2[i];
        distance = std::sqrt(distance);

        const double ekin_v2   = quan.get_kinetic_energy(c2, index_H);
        const double ekin_v1   = quan.get_kinetic_energy(c1, index_H);

        return (ekin_v2 - ekin_v1) / distance;  //Note: This works for electrons and holes
      }

    };

    /** @brief Returns the height of the control box with respect to energy at the provided element (vertex or edge) and energy. */
    template <typename SHEQuantity>
    double box_height(SHEQuantity const & quan,
                      viennagrid_element_id elem,
                      std::size_t index_H)
    {
      //NOTE: Don't try to be clever by being more accurate.
      //      If the box height on the band edge is not taken uniform, carrier conservation is lost!
      double energy_height = 0;

      if (index_H > 0)
        energy_height = std::fabs(quan.get_kinetic_energy(elem, index_H) - quan.get_kinetic_energy(elem, index_H-1));

      if (index_H < quan.get_value_H_size() - 1)
        energy_height += std::fabs(quan.get_kinetic_energy(elem, index_H+1) - quan.get_kinetic_energy(elem, index_H));

      return energy_height / 2.0;
    }


    /** @brief Returns the density of states around a vertex or an edge at total energy specified by index_H. Some averaging is applied near the band edge. */
    template <typename SHEQuantity>
    double averaged_density_of_states(SHEQuantity const & quan,
                                      viennashe::config::dispersion_relation_type const & dispersion,
                                      viennagrid_element_id cell_facet,
                                      std::size_t index_H)
    {

      // Apply some simple averaging over the box center, the upper and the lower energy at the box boundary. Can (should?) be replaced by better integration routines.
      double Z = dispersion.density_of_states(quan.get_kinetic_energy(cell_facet, index_H));
      double num_contributions = 1;

      if (index_H < quan.get_value_H_size() - 1)
      {
        //evaluate density of states at upper boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H + 1)) / 2.0;
        Z += dispersion.density_of_states(kin_energy);
        ++num_contributions;
      }
      if (index_H > 0)
      {
        //evaluate density of states at lower boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H - 1)) / 2.0;
        Z += dispersion.density_of_states(kin_energy);
        ++num_contributions;
      }
      Z /= num_contributions;

      return Z;
    }


    /** @brief Computes \\int v dH using a Simpson rule, where v is velocity
     *
     * Cf. paper by Hong and Jungemann: A fully coupled scheme..., 2009
     *
     * @param dispersion     The dispersion relation used for the SHE simulation
     * @param quan           The SHE quantities
     * @param cell_facet     A topological element which is either a cell or facet
     * @param index_H        The index at which to evaluate the integral
     *
    */
    template <typename SHEQuantity>
    double integral_v(SHEQuantity const & quan,
                      viennashe::config::dispersion_relation_type const & dispersion,
                      viennagrid_element_id cell_facet,
                      std::size_t index_H)
    {

      // Apply some simple averaging over the box center, the upper and the lower energy at the box boundary. Can (should?) be replaced by better integration routines.
      double v = dispersion.velocity(quan.get_kinetic_energy(cell_facet, index_H));
      double num_contributions = 1;

      if (index_H < quan.get_value_H_size() - 1)
      {
        //evaluate density of states at upper boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H + 1)) / 2.0;
        v += dispersion.velocity(kin_energy);
        ++num_contributions;
      }
      if (index_H > 0)
      {
        //evaluate density of states at lower boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H - 1)) / 2.0;
        v += dispersion.velocity(kin_energy);
        ++num_contributions;
      }
      v /= num_contributions;

      return v * box_height(quan, cell_facet, index_H);

    }

    /** @brief Computes \\int v Z dH using a Simpson rule, where v is velocity and Z is the density of states
     *
     * Cf. paper by Hong and Jungemann: A fully coupled scheme..., 2009
     *
     * @param dispersion     The dispersion relation used for the SHE simulation
     * @param quan           The SHE quantities
     * @param cell_facet     A topological element which is either a cell or a facet
     * @param index_H        The index at which to evaluate the integral
     *
    */
    template <typename SHEQuantity>
    double integral_vZ(SHEQuantity const & quan,
                       viennashe::config::dispersion_relation_type const & dispersion,
                       viennagrid_element_id cell_facet,
                       std::size_t index_H)
    {

      // Apply some simple averaging over the box center, the upper and the lower energy at the box boundary. Can (should?) be replaced by better integration routines.
      double kin_energy_center = quan.get_kinetic_energy(cell_facet, index_H);
      double vZ = dispersion.velocity(kin_energy_center) * dispersion.density_of_states(kin_energy_center);
      double num_contributions = 1;

      if (index_H < quan.get_value_H_size() - 1)
      {
        //evaluate density of states at upper boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H + 1)) / 2.0;
        vZ += dispersion.velocity(kin_energy) * dispersion.density_of_states(kin_energy);
        ++num_contributions;
      }
      if (index_H > 0)
      {
        //evaluate density of states at lower boundary of the box:
        double kin_energy = (quan.get_kinetic_energy(cell_facet, index_H)
                            + quan.get_kinetic_energy(cell_facet, index_H - 1)) / 2.0;
        vZ += dispersion.velocity(kin_energy) * dispersion.density_of_states(kin_energy);
        ++num_contributions;
      }
      vZ /= num_contributions;

      return vZ * box_height(quan, cell_facet, index_H);
    }

    /** @brief Returns the averaged kinetic energy in the box centered at a vertex v with total energy index index_H */
    template <typename SHEQuantity>
    double averaged_kinetic_energy(SHEQuantity const & quan,
                                   viennagrid_element_id c,
                                   std::size_t index_H)
    {
      double kin_energy = 0.0;

      kin_energy = quan.get_kinetic_energy(c, index_H);

      if (kin_energy <= 0.0)
      {
        if (index_H < quan.get_value_H_size() - 1)
          return quan.get_kinetic_energy(c, index_H + 1) / 2.0;
      }

      return kin_energy;
    }


  } //namespace she
} //namespace viennashe

#endif

