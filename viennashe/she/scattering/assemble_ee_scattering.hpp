#ifndef VIENNASHE_SHE_ASSEMBLE_EE_SCATTERING_HPP
#define VIENNASHE_SHE_ASSEMBLE_EE_SCATTERING_HPP
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

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/harmonics_coupling.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/filter.hpp"

#include "viennashe/she/assemble_scattering.hpp"
#include "viennashe/she/scattering/all.hpp"
#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/she/postproc/carrier_energy.hpp"
#include "viennashe/she/df_wrappers.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/exception.hpp"


/** @file viennashe/she/scattering/assemble_ee_scattering.hpp
    @brief Implements the assembly for electron-electron scattering
*/

namespace viennashe
{
  namespace she
  {

    template <typename DispersionRelation>
    double ee_scattering_rate(DispersionRelation const & dispersion, double kinetic_energy, double n, double T)
    {
      const double pi   = viennashe::math::constants::pi;
      const double kB   = viennashe::physics::constants::kB;
      const double hbar = viennashe::physics::constants::hbar;
      const double q    = viennashe::physics::constants::q;
      const double eps_Si = viennashe::materials::si::permittivity();
      double lambda_sq = 0.0;
      double scattering_rate = 0.0;

      const double prefactor = (2.0 * pi * n * q * q * q * q) / (hbar * eps_Si * eps_Si);

      double norm_k = dispersion.norm_k(kinetic_energy);

      if (norm_k <= 0.0) //avoid singularity in the following expression
        return 0;

      lambda_sq = (eps_Si * kB * T) / (q * q * n);
      double a = 4.0 * lambda_sq * norm_k * norm_k;

      scattering_rate = prefactor
                        * 0.5 * (std::log(1.0 + a) - (a / (1.0 + a)) ) / 4.0 / std::pow(norm_k, 4.0);

      return scattering_rate;

    }



    /** @brief Interface function for electron-electron scattering.
     *         Differs significantly from ac, op and impurity scattering, thus separate a implementation is used (at least for the moment)
     *
     * @param device  The device description
     * @param conf    The simulator configuration
     * @param matrix  System matrix
     * @param rhs     Load vector
     * @param quan    The current unkown SHE quantity
     * @param quan_old The last computed SHE quantity (probably from the last non-linear iteration step)
    */
    template <typename DeviceType, typename SHEQuantityT, typename MatrixType, typename VectorType>
    void assemble_ee_scattering(DeviceType const & device,
                                viennashe::config const & conf,
                                SHEQuantityT const & quan,
                                SHEQuantityT const & quan_old,
                                MatrixType & matrix, VectorType & rhs
                               )
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

      typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type     CellOnFacetContainer;
      typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                            CellOnFacetIterator;

      typedef viennashe::math::sparse_matrix<double>   CouplingMatrixType;


      // We only assemble EE scattering for electrons
      if ( quan_old.get_carrier_type_id() != viennashe::ELECTRON_TYPE_ID)
        return;


      (void)rhs; // avoid unused parameter warnings

      MeshType const & mesh = device.mesh();

      viennashe::physics::dispersion_proxy dispersion_relation = conf.dispersion_relation(viennashe::ELECTRON_TYPE_ID);

      double damping = 1e-4;  //Empirical fitting factor. Certainly requires more investigations for a general ee-scattering model.

      //
      // Create wrappers for quantities
      //

      viennashe::she::detail::carrier_density_wrapper_by_reference<SHEQuantityT>
          density_n(conf, quan_old);

      viennashe::she::carrier_energy_wrapper<SHEQuantityT>   energy_n(conf, quan_old);

      //
      // Set up scatter matrices
      //
      const std::size_t L_max = static_cast<std::size_t>(conf.max_expansion_order());
      const std::size_t num_harmonics = (L_max+1) * (L_max+1);
      CouplingMatrixType scatter_op_in(num_harmonics, num_harmonics);
      CouplingMatrixType scatter_op_out(num_harmonics, num_harmonics);

      for (std::size_t i = 0; i < (L_max+1) * (L_max+1); ++i)
        scatter_op_out(i,i) += 1.0;
      scatter_op_in(0,0) += 1.0;

      //
      // Step 1: assemble on even nodes:
      //
      //log::debug<log_assemble_ee_scattering>() << "* assembleScatteringOperator(): Assembling on even nodes..." << std::endl;
      CellContainer cells(mesh);
      for (CellIterator cit = cells.begin();
          cit != cells.end();
          ++cit)
      {
        // Step 1: Get carrier density:
        const double carrier_density     = density_n(*cit); // we need the carrier density in a single valley
        const double kinetic_energy_star = energy_n(*cit);  // mean particle energy

        //
        // Iterate over relevant energies (at index_H == 0 there is no unknown...)
        for (std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H)
        {

          const long row_index = quan.get_unknown_index(*cit, index_H);

          if (row_index < 0)
            continue;

          const double kinetic_energy = quan.get_kinetic_energy(*cit, index_H);

          if (kinetic_energy <= 0.0)
            continue;

          if (carrier_density <= 0)
          {
            log::warn() << "* assemble_ee_scattering(): Warning: Carrier_density zero!" << std::endl;
          }
          if (kinetic_energy_star <= 0)
          {
            log::warn() << "* assemble_ee_scattering(): Warning:  kinetic_energy_star zero or less ! kinetic_energy_star = " << kinetic_energy_star << std::endl;
          }

          const double box_volume        = viennagrid::volume(*cit);
          const long expansion_order_mid = static_cast<long>(quan.get_expansion_order(*cit, index_H));
          const double energy_height     = box_height(quan, *cit, index_H);


          //////////////////////////////////////////////////////////////////////////////////////////////////
          // Step 2: Determine E* and f(E*):
          std::size_t Estar_index = 1;
          for ( ; Estar_index < quan_old.get_value_H_size() - 1; ++Estar_index)
          {
            if (quan_old.get_kinetic_energy(*cit, Estar_index) > kinetic_energy_star)
              break;
          }

          double f_Estar = quan_old.get_values(*cit, Estar_index)[0] / carrier_density;  //normalize f such that \int f dk^3 = 1

          { // account for DOS ... we now assemble for g = f*Z
            const double dos_norm = dispersion_relation.density_of_states(quan_old.get_kinetic_energy(*cit, Estar_index));
            if (dos_norm <= 0.0)
            {
              f_Estar = 0;
              std::cout << "kinetic_energy = " << quan_old.get_kinetic_energy(*cit, Estar_index) << std::endl;
              std::cout << "dos_norm       = " << dos_norm << std::endl;
            }
            else f_Estar = f_Estar / dos_norm; // account for DOS ... we now assemble for g = f*Z
          }

          if (f_Estar <= 0)
          {
            log::warn() << "* assemble_ee_scattering(): Warning:  f_Estar is smaller or equal to zero!" << std::endl;
          }

          // Step 3: Run an integration over energy to obtain scattering contributions:
          double out_scattering_coeff = 0.0;
          double kinetic_energy_max = quan.get_kinetic_energy(*cit, quan_old.get_value_H_size() - 2);
          viennashe::math::SphericalHarmonic Y_00(0, 0);
          double prefactor = pow(1.0 / Y_00(0, 0), 3) / 2.0; //prefactor for the integrals

          for (std::size_t index_H_prime = 1; index_H_prime < quan.get_value_H_size() - 1; ++index_H_prime)
          {
            const double kinetic_energy_prime = quan.get_kinetic_energy(*cit, index_H_prime);

            if (kinetic_energy_prime <= 0)
              continue;

            //store E + E* - E' in a separte variable:
            const double energy_argument = kinetic_energy + kinetic_energy_star - kinetic_energy_prime;
            if (energy_argument < 0) //since E' increases monotonically, energy_argument will only get smaller -> break for-loop
              break;

            if (energy_argument > kinetic_energy_max)
              continue;

            const double c_ee = ee_scattering_rate(dispersion_relation,
                                                   std::abs(kinetic_energy - kinetic_energy_prime),
                                                   carrier_density,
                                                   device.get_lattice_temperature(*cit));

            //
            // In-scattering part: Energies couple
            //
            double f_00_prime = (quan_old.get_values(*cit, index_H_prime)[0]) / carrier_density;  //normalization to \int f dk^3 = 1
            { // account for DOS ... we now assemble for g = f*Z
              const double dos_norm = dispersion_relation.density_of_states(kinetic_energy_prime);
              if (dos_norm <= 0.0) f_00_prime = 0;
              else f_00_prime = f_00_prime / dos_norm;
            }

            // find unknown index for g_00(E + E* - E'):
            std::size_t energy_index = 1;
            for ( ; energy_index < quan.get_value_H_size() - 1; ++energy_index)
            {
              if (quan.get_kinetic_energy(*cit, energy_index) > energy_argument)
                break;
            }

            if (energy_index == quan.get_value_H_size() - 2) //E + E* - E' is outside the simulation mesh
              continue;

            double in_scattering_coeff = prefactor
                                        * c_ee
                                        * dispersion_relation.density_of_states(kinetic_energy_prime)
                                        * dispersion_relation.density_of_states(energy_argument)
                                        * energy_height; // \int_E h(E) dE   is replaced by    \sum_Ei h(Ei) Delta(E_i)

            // write to system matrix: (mind negative sign of in-scattering)
            if (quan.get_unknown_index(*cit, energy_index) >= 0)
              viennashe::util::add_block_matrix(matrix,
                                                std::size_t(row_index), std::size_t(quan.get_unknown_index(*cit, energy_index)),
                                                - damping * in_scattering_coeff * f_00_prime  * box_volume * energy_height,
                                                scatter_op_in,
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                                viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::EVEN_HARMONICS_ITERATION_ID)
                                              );

            //
            // Out-scattering part:
            //
            out_scattering_coeff += prefactor
                                    * c_ee
                                    * dispersion_relation.density_of_states(kinetic_energy_prime)
                                    * dispersion_relation.density_of_states(energy_argument)
                                    * energy_height; // \int_E h(E) dE   is replaced by    \sum_Ei h(Ei) Delta(E_i)
          }

          //write out-scattering terms:
          double coeff = damping * out_scattering_coeff * f_Estar  * box_volume * energy_height;
          if (coeff > 0)
            viennashe::util::add_block_matrix(matrix,
                                              std::size_t(row_index), std::size_t(row_index),
                                              coeff,
                                              scatter_op_out,
                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::EVEN_HARMONICS_ITERATION_ID),
                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::EVEN_HARMONICS_ITERATION_ID)
                                            );
        } //index_H
      } //for vertices


      //
      // Step 2: assemble on odd nodes:
      //
      //log::debug<log_assemble_ee_scattering>() << "* assembleScatteringOperator(): Assembling on odd nodes..." << std::endl;
      FacetContainer facets(mesh);
      for (FacetIterator fit = facets.begin();
           fit != facets.end();
           ++fit)
      {
        CellOnFacetContainer cells_on_facet(device.mesh(), fit.handle());

        if (cells_on_facet.size() < 2)
          continue;

        CellOnFacetIterator cofit = cells_on_facet.begin();
        CellType const & c1 = *cofit;
        ++cofit;
        CellType const & c2 = *cofit;

        double connection_len = viennagrid::norm_2(viennagrid::centroid(c1) - viennagrid::centroid(c2));

        //
        // Iterate over relevant energies (at index_H == 0 there is no unknown...)
        for (std::size_t index_H = 1; index_H < quan.get_value_H_size() - 1; ++index_H)
        {
          const long row_index = quan.get_unknown_index(*fit, index_H);
          if (row_index < 0)
            continue;


          const double box_volume = viennagrid::volume(*fit) * connection_len;

          const long expansion_order_mid = static_cast<long>(quan.get_expansion_order(*fit, index_H));

          //
          // Compute density of states at initial and final states
          //
          const double kinetic_energy_at_edge = quan.get_kinetic_energy(*fit, index_H);
          const double energy_height          = box_height(quan, *fit, index_H);

          //////////////////////////////////////////////////////////////////////////////////////////////////

          //std::cout << density_n(v1) << " * " <<  density_n(v2) << std::endl;
          // Step 1: Get carrier density:
          double carrier_density = std::sqrt(density_n(c1) * density_n(c2));
          if (carrier_density <= 0)
            log::warn() << "* assemble_ee_scattering(): Warning: carrier_density zero!" << std::endl;

          //std::cout << energy_n(v1) << " * " <<  energy_n(v2) << std::endl;
          // Step 2: Determine E* and f(E*):
          double kinetic_energy_star = (energy_n(c1) + energy_n(c2)) / 2.0; //mean particle energy
          if (kinetic_energy_star <= 0)
            log::warn() << "* assemble_ee_scattering(): Warning: kinetic_energy_star zero!" << std::endl;

          std::size_t Estar_index = 1;
          for ( ; Estar_index < quan_old.get_value_H_size(); ++Estar_index)
          {
            if (quan_old.get_kinetic_energy(*fit, Estar_index) > kinetic_energy_star)
              break;
          }



          double f_Estar = 0;

          if ( quan_old.get_unknown_mask(c1) ) f_Estar += quan_old.get_boundary_value(c1, Estar_index);
          else f_Estar += quan_old.get_values(c1, Estar_index)[0];

          if ( quan_old.get_unknown_mask(c2) ) f_Estar += quan_old.get_boundary_value(c2, Estar_index);
          else f_Estar += quan_old.get_values(c2, Estar_index)[0];

          f_Estar = f_Estar / (2.0 * carrier_density);

          //f_Estar = (  quan_old.get_values(v1, Estar_index)[0]
          //                  + quan_old.get_values(v2, Estar_index)[0] ) / (2.0 * carrier_density);  //normalize f such that \int f dk^3 = 1

          {
            const double dos_norm = dispersion_relation.density_of_states(quan_old.get_kinetic_energy(*fit, Estar_index));
            if (dos_norm <= 0.0) f_Estar = 0.0;
            else f_Estar = f_Estar / dos_norm;
          }

          // Step 3: Run an integration over energy to obtain scattering contributions:
          double kinetic_energy_max = quan.get_kinetic_energy(*fit, quan_old.get_value_H_size() - 2);
          double out_scattering_coeff = 0.0;
          viennashe::math::SphericalHarmonic Y_00(0, 0);
          double prefactor = pow(1.0 / Y_00(0, 0), 3); //prefactor for the integrals

          for (std::size_t index_H_prime = 1; index_H_prime < quan.get_value_H_size() - 1; ++index_H_prime)
          {
            double kinetic_energy_prime = quan.get_kinetic_energy(*fit, index_H_prime);
            if (kinetic_energy_prime < 0)
              continue;

            //store E + E* - E' in a separte variable:
            const double energy_argument = kinetic_energy_at_edge + kinetic_energy_star - kinetic_energy_prime;
            if (energy_argument < 0) //since E' increases monotonically, energy_argument will only get smaller -> break for-loop
              break;

            if (energy_argument > kinetic_energy_max)
              continue;

            const double c_ee = ee_scattering_rate( dispersion_relation,
                                                    std::abs(kinetic_energy_at_edge - kinetic_energy_prime),
                                                    carrier_density,
                                                    (device.get_lattice_temperature(c1) + device.get_lattice_temperature(c2)) / 2.0
                                                  );

            //
            // In-scattering part: Nothing to be done, because couples only f_{0,0}
            //

            //
            // Out-scattering part:
            //
            out_scattering_coeff += prefactor
                                    * c_ee
                                    * dispersion_relation.density_of_states(kinetic_energy_prime)
                                    * dispersion_relation.density_of_states(energy_argument)
                                    * energy_height; // \int_E h(E) dE   is replaced by    \sum_Ei h(Ei) Delta(E_i)
          }

          //write out-scattering terms:
          const double coeff = damping * out_scattering_coeff * f_Estar * box_volume * energy_height;
          if (coeff > 0)
            viennashe::util::add_block_matrix(matrix,
                                              std::size_t(row_index), std::size_t(row_index),
                                              coeff,
                                              scatter_op_out,
                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::ODD_HARMONICS_ITERATION_ID),
                                              viennashe::math::spherical_harmonics_iterator(expansion_order_mid, viennashe::math::ODD_HARMONICS_ITERATION_ID)
                                            );

        } //index_H
      } //for edges



    } // assemble_ee

  } //namespace she
} //namespace viennashe

#endif

