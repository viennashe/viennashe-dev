#ifndef VIENNASHE_DD_ASSEMBLE_HPP
#define VIENNASHE_DD_ASSEMBLE_HPP

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
#include <fstream>
#include <vector>

// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/mesh/coboundary_iteration.hpp"
#include "viennagrid/algorithm/volume.hpp"

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/device.hpp"
#include "viennashe/math/bernoulli.hpp"
#include "viennashe/math/linalg_util.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/util/dual_box_flux.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/log_keys.h"

#include "viennashe/scharfetter_gummel.hpp"

#include "viennashe/accessors.hpp"
#include "viennashe/she/timestep_quantities.hpp"
#include "viennashe/she/assemble_all.hpp"

#include "viennashe/phonon/joule_heating.hpp"


/** @file viennashe/assemble.hpp
    @brief Dimension-independent assembling routines for the Drift-diffusion model using a Newton scheme or Gummel scheme for self-consistency
*/

namespace viennashe
{

  namespace detail
  {
    /**
     * @brief Assembles, per cell, the carrier concentration coupling on the RHS of Poisson's equation
     * @param conf The simulator configuration
     * @param cell The cell on which to assemble
     * @param np_density The carrier density
     * @param f_np The SHE guess/solution
     * @param A The system matrix
     * @param row_index The row index associated with the potential unkown and the vertex
     * @param box_volume The volume of the current finite box associated with the vertex
     *
     */
    template <typename CellType, typename SpatialUnknownT, typename SHEUnknownT, typename MatrixType>
    void assemble_poisson_carrier_coupling(viennashe::config const & conf,
                                           CellType const & cell,
                                           SpatialUnknownT const & np_density,
                                           SHEUnknownT const & f_np,
                                           MatrixType & A, std::size_t row_index, double box_volume)
    {
      const double polarity = (f_np.get_carrier_type_id() == ELECTRON_TYPE_ID) ? -1.0 : 1.0;
      equation_id equ_id = (f_np.get_carrier_type_id() == ELECTRON_TYPE_ID) ? conf.get_electron_equation() : conf.get_hole_equation();

      if (equ_id == EQUATION_CONTINUITY)  //Drift-Diffusion: Coupling term -|q|*n for electrons, |q|*p for holes
      {
        const long np_index = np_density.get_unknown_index(cell);  // This is a hack and needs to be resolved with ViennaMath!
        if (np_index >= 0)  //here is a Dirichlet boundary condition or no simulation region
          A(row_index, std::size_t(np_index)) = polarity * viennashe::physics::constants::q * box_volume;
      }
      else //SHE: -|q| * integral_{energies} f_n(E) * Z(E) dE for electrons, |q| * integral_{energies} f_p(E) * Z(E) dE for holes
      {
        // matrix entries (note: the following piece of code pretty much resembles the one in postproc/carrier_density.hpp - which is required to be consistent)
        for (std::size_t index_H = 1; index_H < f_np.get_value_H_size() - 1; ++index_H)
        {
          long col_index = f_np.get_unknown_index(cell, index_H);
          if (col_index >= 0)
          {
            double energy_mid = f_np.get_kinetic_energy(cell, index_H);
            double energy_upper =                   (f_np.get_kinetic_energy(cell, index_H + 1) + energy_mid) / 2.0;
            double Z            = viennashe::she::averaged_density_of_states(f_np, conf.dispersion_relation(f_np.get_carrier_type_id()), cell, index_H);

            if (energy_upper <= 0.0) //nothing to do if in band gap
              continue;

            double box_height = viennashe::she::box_height(f_np, cell, index_H);

            // Check Y_{0,0} !!
            A(row_index, std::size_t(col_index)) =  polarity
                                                     * viennashe::physics::constants::q
                                                     * box_volume
                                                     * Z * conf.dispersion_relation(f_np.get_carrier_type_id()).symmetry_factor() * box_height;  // Check Y_{0,0} !!
          }
        } // for index_H
      }
    }

    /**
     * @brief Returns the requested carrier concentration, calculated from a SHE or DD solution, for the RHS of Poisson's equation
     * @param conf The simulator configuration
     * @param cell The vertex on which the concentration is requested
     * @param np_density The spatial unkown for n or p (dummy if not used)
     * @param f_np The SHE solution (dummy if not used)
     * @return The requested carrier concentration
     */
    template <typename DeviceType, typename CellType, typename SpatialUnknownT, typename SHEUnknownT>
    double get_carrier_density_for_poisson(viennashe::config const & conf,
                                           CellType          const & cell,
                                           SpatialUnknownT   const & np_density,
                                           SHEUnknownT       const & f_np)
    {
      const equation_id equ_id    = (f_np.get_carrier_type_id() == ELECTRON_TYPE_ID) ? conf.get_electron_equation() : conf.get_hole_equation();
      const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

      if (equ_id == EQUATION_CONTINUITY || !with_full_newton)
      {
        return np_density.get_value(cell);
      }
      else
      {
        viennashe::she::carrier_density_wrapper<SHEUnknownT> density_wrapper(conf, f_np);
        return density_wrapper(cell);
      }
    }

  }

  /** @brief Assembles the potential block (i.e. the linear equations obtained from the discretization of the Poisson equation) in the Jacobian matrix.
   *
   *  Note that this routine is used both for a Newton method and for a Gummel/Picard-iteration method for the solution of the nonlinear system of drift-diffusion equations.
   *
   * @tparam DeviceType The ViennaSHE device type
   * @tparam MatrixType The matrix type
   * @tparam VectorType The vector type
   *
   * @param device      The device to be simulated
   * @param quantities  Quantity descriptors and accessors
   * @param conf        The simulator configuration
   * @param A           The system matrix
   * @param b           The right hand side
   */
  template <typename DeviceType,
            typename MatrixType,
            typename VectorType>
  void assemble_poisson(DeviceType const & device,
                        viennashe::she::timestep_quantities<DeviceType> const & quantities,
                        viennashe::config const & conf,
                        MatrixType & A,
                        VectorType & b)
  {
    typedef typename DeviceType::mesh_type           MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
    typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

    typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_quantity_type      SpatialUnknownType;
    typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_she_quantity_type  SHEUnknownType;

    //
    // Poisson equation:  + eps * laplace psi - |q| * (n - p - doping) = 0
    //   Linearisation (Gummel):  + eps * laplace dpsi - q * (n+p)/VT * dpsi = - eps * laplace psi + q * (n - p - doping)
    // Note that the term q * (n+p)/VT is for Gummel damping
    //

    MeshType const & mesh = device.mesh();

    bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

    // Set up filters:
    permittivity_accessor<DeviceType>              permittivity(device);
    fixed_charge_accessor<DeviceType>              fixed_charge(device);

    // Quantity lookup:
    SpatialUnknownType const & potential = quantities.get_unknown_quantity(viennashe::quantity::potential());
    SpatialUnknownType const & n_density = quantities.get_unknown_quantity(viennashe::quantity::electron_density());
    SpatialUnknownType const & p_density = quantities.get_unknown_quantity(viennashe::quantity::hole_density());
    SHEUnknownType     const & f_n       = quantities.electron_distribution_function();
    SHEUnknownType     const & f_p       = quantities.hole_distribution_function();

    CellContainer cells(mesh);
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {

      const long row_index2 = potential.get_unknown_index(*cit);
      if (row_index2 < 0)  //here is a Dirichlet boundary condition
        continue;
      const std::size_t row_index = std::size_t(row_index2);

      PointType centroid_cell = viennagrid::centroid(*cit);

      const double potential_center = potential.get_value(*cit);
      const double permittivity_center = permittivity(*cit);


      b[row_index] = 0;
      A(row_index, row_index) = 0;

      //
      //   eps * laplace psi
      //
      FacetOnCellContainer facets(*cit);
      for (FacetOnCellIterator focit = facets.begin();
          focit != facets.end();
          ++focit)
      {
        CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *focit, *cit);

        if (!other_cell_ptr) continue;  //Facet is on the boundary of the simulation domain -> homogeneous Neumann conditions

        PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
        PointType cell_connection = centroid_other_cell - centroid_cell;
        PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
        PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);

        const long col_index        = potential.get_unknown_index(*other_cell_ptr);
        const double connection_len = viennagrid::norm_2(cell_connection);
        const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);
        const double potential_outer = potential.get_value(*other_cell_ptr);

        // off-diagonal contribution
        if (col_index >= 0)
        {
          PointType facet_center = viennagrid::centroid(*focit);  //TODO: Use intersection of facet plane with connection
          double connection_in_cell = viennagrid::norm(facet_center - centroid_cell);
          double connection_in_other_cell = viennagrid::norm(facet_center - centroid_other_cell);
          const double permittivity_mean = (connection_in_cell + connection_in_other_cell) /
                                           (connection_in_cell/permittivity_center + connection_in_other_cell/permittivity(*other_cell_ptr));

          A(row_index, std::size_t(col_index))  = permittivity_mean * weighted_interface_area / connection_len;
          A(row_index, row_index) -= permittivity_mean * weighted_interface_area / connection_len;
          b[row_index]            -= permittivity_mean * weighted_interface_area * (potential_outer - potential_center) / connection_len;
        }
        else if (potential.get_boundary_type(*other_cell_ptr) == BOUNDARY_DIRICHLET)
        {
          A(row_index, row_index) -= permittivity_center * weighted_interface_area / connection_len;
          b[row_index]            -= permittivity_center * weighted_interface_area / connection_len * (potential.get_boundary_value(*other_cell_ptr) - potential_center);
        }
      } //for facets

      const double cell_volume = viennagrid::volume(*cit);

      const double value_n = detail::get_carrier_density_for_poisson<DeviceType>(conf, *cit, n_density, f_n);
      const double value_p = detail::get_carrier_density_for_poisson<DeviceType>(conf, *cit, p_density, f_p);

      //
      // Gummel-Damping:  extra volume term q * (n + p) / V_T
      //
      if (!with_full_newton && viennashe::materials::is_semiconductor(device.get_material(*cit)))
      {
        const double VT = viennashe::physics::get_thermal_potential(device.get_lattice_temperature(*cit));

        A(row_index, row_index) -= cell_volume * viennashe::physics::constants::q * (value_n + value_p) / VT;
      }

      //
      // - |q| * (n - p - doping)
      //

      // carrier contribution to Jacobian matrix  - |q| * n + |q| * p
      if (with_full_newton)
      {
        if (conf.with_electrons())
          detail::assemble_poisson_carrier_coupling(conf, *cit, n_density, f_n, A, row_index, cell_volume);
        if (conf.with_holes())
          detail::assemble_poisson_carrier_coupling(conf, *cit, p_density, f_p, A, row_index, cell_volume);
      }


      // residual contribution
      const double doping_n = device.get_doping_n(*cit);
      const double doping_p = device.get_doping_p(*cit);

      if (viennashe::materials::is_semiconductor(device.get_material(*cit)))
      {
        b[row_index] += cell_volume * viennashe::physics::constants::q * (value_n - doping_n - value_p + doping_p);
        // *Fixed* charges
        b[row_index] +=  fixed_charge(*cit) * cell_volume; // Charge has to be in As
      }

      if (conf.with_traps() && conf.with_trap_selfconsistency()
           && viennashe::materials::is_semiconductor(device.get_material(*cit))) // *Trapped* charges
      {
        typedef typename DeviceType::trap_level_container_type     trap_level_container_type;
        typedef typename trap_level_container_type::const_iterator trap_iterator_type;

        trap_level_container_type const & traps = device.get_trap_levels(*cit);

        const std::size_t num_trap_unknowns = quantities.num_trap_unknown_indices(*cit);

        if (num_trap_unknowns != traps.size())
          throw viennashe::invalid_value_exception("The number of traps configured in the device does not match the number of unknowns for traps!", static_cast<double>(num_trap_unknowns));

        if (num_trap_unknowns != 0)
        {
          std::size_t index = 0;
          for (trap_iterator_type tit = traps.begin(); tit != traps.end(); ++tit, ++index)
          {
            const double occupancy = quantities.trap_occupancy(*cit, index);
            // Charge has to be in As => sign * |q| * f_T * N_T * V_i
            b[row_index] +=  tit->charge_sign() * viennashe::physics::constants::q * occupancy * tit->density() * cell_volume ;
          } // for trap levels
        }

      } // trapped charges


    } //for cells

  } //assemble_poisson



  /**
   * Assembles the Drift Diffusion (DD) equation set (div J = +-q R) for the given carrier type
   * @param device The device
   * @param quantities The unkown quantities
   * @param conf The simulator configuration
   * @param ctype The carrier type (normally electrons or holes)
   * @param A The system matrix
   * @param b The right hand side (RHS)
   */
  template <typename DeviceType,
            typename MatrixType,
            typename VectorType>
  void assemble_dd(DeviceType const & device,
                   viennashe::she::timestep_quantities<DeviceType> const & quantities,
                   viennashe::config const & conf,
                   carrier_type_id ctype,
                   MatrixType & A,
                   VectorType & b)
  {
    typedef typename DeviceType::mesh_type           MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
    typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

    typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_quantity_type      SpatialUnknownType;
    typedef typename viennashe::contact_carrier_density_accessor<DeviceType> bnd_carrier_accessor;

    MeshType const & mesh = device.mesh();

    bnd_carrier_accessor bnd_carrier_density(device, ctype);

    const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

    // Quantity lookup:
    SpatialUnknownType const & potential       = quantities.get_unknown_quantity(viennashe::quantity::potential());
    SpatialUnknownType const & carrier_density = (ctype == ELECTRON_TYPE_ID) ? quantities.get_unknown_quantity(viennashe::quantity::electron_density())
                                                                             : quantities.get_unknown_quantity(viennashe::quantity::hole_density());
    SpatialUnknownType const & quantum_corr    = (ctype == ELECTRON_TYPE_ID) ? quantities.get_unknown_quantity(viennashe::quantity::density_gradient_electron_correction())
                                                                             : quantities.get_unknown_quantity(viennashe::quantity::density_gradient_hole_correction());

    double mobility = 1.0; //[KR] TODO: Add better mobility model here. Doesn't really matter as long as there is no recombination...

    scharfetter_gummel flux_approximator(ctype);
    scharfetter_gummel_dVi flux_approximator_dVi(ctype);
    scharfetter_gummel_dVj flux_approximator_dVj(ctype);

    CellContainer cells(mesh);
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {
      //log::debug<log_continuity_solver>() << "* electron_solver::assemble(): Iterating over vertex: " << vit->id() << std::endl;

      const long row_index2 = carrier_density.get_unknown_index(*cit);
      if (row_index2 < 0)
        continue;
      const std::size_t row_index = std::size_t(row_index2);

      //get potential in box center:
      const long pot_index_center  = potential.get_unknown_index(*cit);
      double potential_center      = potential.get_value(*cit);

      if (conf.with_quantum_correction())
        potential_center += quantum_corr.get_value(*cit);

      const double carrier_center = carrier_density.get_value(*cit);

      PointType centroid_cell = viennagrid::centroid(*cit);

      A(row_index, row_index) = 0.0; //make sure that there is no bogus in the diagonal
      b[row_index]            = 0.0;

      FacetOnCellContainer facets(*cit);
      for (FacetOnCellIterator focit = facets.begin();
          focit != facets.end();
          ++focit)
      {
        CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *focit, *cit);

        if (!other_cell_ptr) continue;

        const double T = 0.5 * (device.get_lattice_temperature(*cit) + device.get_lattice_temperature(*other_cell_ptr));

        if ( (carrier_density.get_unknown_mask(*other_cell_ptr) || carrier_density.get_boundary_type(*other_cell_ptr) == BOUNDARY_DIRICHLET) )
        {
          PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
          PointType cell_connection = centroid_other_cell - centroid_cell;
          PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
          PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);

          const double connection_len = viennagrid::norm_2(centroid_cell - centroid_other_cell);
          const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);
          double potential_outer = potential.get_value(*other_cell_ptr);
          if (conf.with_quantum_correction())
            potential_outer += quantum_corr.get_value(*other_cell_ptr);

          const long pot_index_other  = potential.get_unknown_index(*other_cell_ptr);

          //
          // div() part of continuity equation:
          //
          const long np_col_index = carrier_density.get_unknown_index(*other_cell_ptr);
          if (np_col_index > -1) // no bc
          {
            A(row_index, std::size_t(np_col_index)) = weighted_interface_area * flux_approximator(0.0, 1.0,
                                                                                                  potential_center, potential_outer,
                                                                                                  connection_len, mobility, T);
          }

          A(row_index, row_index) += weighted_interface_area * flux_approximator(1.0, 0.0,
                                                                                 potential_center, potential_outer,
                                                                                 connection_len, mobility, T);

          //
          // Residual contribution
          //
          const double carrier_outer = viennashe::materials::is_conductor(device.get_material(*other_cell_ptr))
                                       ? bnd_carrier_density(*cit, device.get_lattice_temperature(*other_cell_ptr))
                                       : carrier_density.get_value(*other_cell_ptr);

          b[row_index] -= weighted_interface_area * flux_approximator(carrier_center, carrier_outer,
                                                                      potential_center, potential_outer,
                                                                      connection_len, mobility, T);

          if (!with_full_newton)  //Gummel-type iteration ignores cross-couplings
            continue;

          //
          // NEWTON-SOLVER: couplings for Jacobi matrix
          //

          // DG contributions
          if (conf.with_quantum_correction())
          {
            const long corrpot_center_index = quantum_corr.get_unknown_index(*cit);
            const long corrpot_outer_index  = quantum_corr.get_unknown_index(*other_cell_ptr);
            if (corrpot_center_index >= 0)
            {
              const double gamma_outer  = quantum_corr.get_value(*other_cell_ptr);
              const double gamma_center = quantum_corr.get_value(*cit);
              A(row_index, std::size_t(corrpot_center_index)) = weighted_interface_area *
                                                                flux_approximator_dVj(carrier_center, carrier_outer,
                                                                                      gamma_center, gamma_outer,
                                                                                      connection_len, mobility, T);
            }
            if (corrpot_outer_index >= 0)
            {
              const double gamma_outer  = quantum_corr.get_value(*other_cell_ptr);
              const double gamma_center = quantum_corr.get_value(*cit);
              A(row_index, std::size_t(corrpot_center_index)) = weighted_interface_area *
                                                                flux_approximator_dVj(carrier_center, carrier_outer,
                                                                                      gamma_center, gamma_outer,
                                                                                      connection_len, mobility, T);
            }
          }

          // off-diagonal contribution
          if (pot_index_other >= 0)
          {
            A(row_index, std::size_t(pot_index_other)) = weighted_interface_area *
                                                         flux_approximator_dVj(carrier_center, carrier_outer,
                                                                               potential_center, potential_outer,
                                                                               connection_len, mobility, T);
          }

          // 'diagonal' contribution
          A(row_index, std::size_t(pot_index_center)) += weighted_interface_area *
                                                         flux_approximator_dVi(carrier_center, carrier_outer,
                                                                               potential_center, potential_outer,
                                                                               connection_len, mobility, T);
        }

      } //for edges

    } //for vertices
  }


  namespace detail
  {
    /**
     * @brief Returns the quantum correction potential flux between two points (Density Gradient)
     * @param E The band edge energy
     * @param pot_i Potential at the source
     * @param gamma_i Correction potential at the source
     * @param T_i Temperature (lattice) at the source
     * @param pot_j Potential at the sink
     * @param gamma_j Correction potential at the sink
     * @param T_j Temperature at the sink
     * @return The flux
     */
    inline double density_gradient_flux(double E, double pot_i, double gamma_i, double T_i, double pot_j, double gamma_j, double T_j)
    {
      const double kB = viennashe::physics::constants::kB;
      const double q = viennashe::physics::constants::q;
      const double kbT_i = (kB * T_i) * 2.0;
      const double kbT_j = (kB * T_j) * 2.0;
      return ( (q * (pot_j + gamma_j) - E) / kbT_j) - ((q * (pot_i + gamma_i) - E) / kbT_i);
    }
  }

  /**
   * @brief Assembles the density gradient equation set for the given carrier type
   * @param device The device
   * @param quantities The unkown quantities and required known ones
   * @param conf The simulator configruation
   * @param ctype The carrier type for which to assemble
   * @param A The system matrix
   * @param b The right hand side (RHS)
   *
   */
  template <typename DeviceType, typename MatrixType, typename VectorType>
  void assemble_density_gradient(DeviceType const & device,
                                 viennashe::she::timestep_quantities<DeviceType> const & quantities,
                                 viennashe::config const & conf,
                                 carrier_type_id ctype,
                                 MatrixType & A,
                                 VectorType & b)
  {
    typedef typename DeviceType::mesh_type           MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
    typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

    typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_quantity_type      SpatialUnknownType;

    MeshType const & mesh = device.mesh();

    const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

    // Quantity lookup:
    SpatialUnknownType const & potential       = quantities.get_unknown_quantity(viennashe::quantity::potential());
    SpatialUnknownType const & quantum_corr    = (ctype == ELECTRON_TYPE_ID) ? quantities.get_unknown_quantity(viennashe::quantity::density_gradient_electron_correction())
                                                                             : quantities.get_unknown_quantity(viennashe::quantity::density_gradient_hole_correction());

    const double kB     = viennashe::physics::constants::kB;
    const double q      = viennashe::physics::constants::q;
    const double lambda = conf.density_gradient(ctype).lambda();

    const double m0       = viennashe::physics::constants::mass_electron;
    const double coeff_b  = (viennashe::physics::constants::hbar * viennashe::physics::constants::hbar) / (6.0 * lambda * m0 * q);

    const double nominal_band_edge = viennashe::physics::get_band_edge(ctype);

    CellContainer cells(mesh);
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {
      const long row_index2 = quantum_corr.get_unknown_index(*cit);
      if ( row_index2 < 0 ) //here is a Dirichlet boundary condition
        continue;
      const std::size_t row_index = std::size_t(row_index2);

      const double potential_center = potential.get_value(*cit);
      const double gamma_center = quantum_corr.get_value(*cit);

      // Reset matrix and rhs for the current index
      A(row_index, row_index) = 0.0;
      b[row_index] = 0.0;

      PointType centroid_cell = viennagrid::centroid(*cit);
      const double box_volume = viennagrid::volume(*cit);

      FacetOnCellContainer facets(*cit);
      for (FacetOnCellIterator focit = facets.begin();
          focit != facets.end();
          ++focit)
      {
        CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *focit, *cit);

        if (!other_cell_ptr) continue;

        const double T = 0.5 * (device.get_lattice_temperature(*cit) + device.get_lattice_temperature(*other_cell_ptr));

        PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
        PointType cell_connection = centroid_other_cell - centroid_cell;
        PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
        PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);

        const double connection_len = viennagrid::norm_2(centroid_cell - centroid_other_cell);
        const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);

        const long  col_index = quantum_corr.get_unknown_index(*other_cell_ptr);

        const double c_ij            = coeff_b * (weighted_interface_area / connection_len);
        const double potential_outer = potential.get_value(*other_cell_ptr);

        // off-diagonal contribution
        if ( col_index >= 0 )
        {
          const double gamma_outer = quantum_corr.get_value(*other_cell_ptr);

          A(row_index, std::size_t(col_index)) -= c_ij * q / (kB * T * 2.0);
          A(row_index, row_index) += c_ij * q / (kB * T * 2.0);

          b[row_index] += c_ij * detail::density_gradient_flux(nominal_band_edge, potential_center, gamma_center, T, potential_outer, gamma_outer, T);
        }
        else if ( quantum_corr.get_boundary_type(*other_cell_ptr) == BOUNDARY_DIRICHLET)
        {
          const double gamma_outer = quantum_corr.get_boundary_value(*other_cell_ptr);

          A(row_index, row_index) += c_ij * q / (kB * T * 2.0);

          b[row_index] += c_ij * detail::density_gradient_flux(nominal_band_edge, potential_center, gamma_center, T, potential_outer, gamma_outer, T);
        }

        //
        // CHECK THIS ...
        if ( with_full_newton ) //Gummel-type iteration ignores cross-couplings
        {
          const long pot_index_outer  = potential.get_unknown_index(*other_cell_ptr);
          const long pot_index_center = potential.get_unknown_index(*cit);

          // off-diagonal contribution
          if ( pot_index_outer >= 0 )
            A(row_index, std::size_t(pot_index_outer)) -= c_ij * q / (kB * T * 2.0);

          // 'diagonal' contribution
          A(row_index, std::size_t(pot_index_center)) += c_ij * q / (kB * T * 2.0);
        }

        if ( quantum_corr.get_boundary_type(*other_cell_ptr) == BOUNDARY_ROBIN )
        {
          // Get boundary coefficients
          robin_boundary_coefficients<double> robin_coeffs = quantum_corr.get_boundary_values(*other_cell_ptr);
          const double alpha = robin_coeffs.alpha;
          const double beta  = robin_coeffs.beta;

          // Apply Robin-condition
          A(row_index, row_index) -= alpha * weighted_interface_area ;
          b[row_index]            += (beta + alpha * gamma_center) * weighted_interface_area;
        }
      } //for facets

      // volume contributions
      A(row_index, row_index) += box_volume;
      b[row_index] -= box_volume * gamma_center;

    } //for vertices

  }

  /**
   * @brief Assembles the heat diffusion equation (HDE)
   * @param device The device
   * @param quantities The unknown (lattice temperature) and required known quantities
   * @param conf The simulator configuration
   * @param A The system matrix
   * @param b The right hand side
   */
  template <typename DeviceType, typename MatrixType, typename VectorType>
  void assemble_heat(DeviceType const & device,
                     viennashe::she::timestep_quantities<DeviceType> const & quantities,
                     viennashe::config const & conf,
                     MatrixType & A,
                     VectorType & b)
  {
    typedef typename viennashe::she::timestep_quantities<DeviceType> QuantitiesType;
    typedef typename DeviceType::mesh_type           MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
    typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

    typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

    typedef typename viennashe::she::timestep_quantities<DeviceType>::unknown_quantity_type      SpatialUnknownType;

    //
    // Poisson equation for heat:  + kappa * laplace T = Q
    //

    MeshType const & mesh = device.mesh();

    diffusivity_accessor<DeviceType>  diffusivity(device);
    viennashe::hde::power_density_accessor<DeviceType, QuantitiesType> quan_power_density(device, quantities, conf);

    // Quantity lookup:
    SpatialUnknownType const & lattice_temperature = quantities.get_unknown_quantity(viennashe::quantity::lattice_temperature());

    CellContainer cells(mesh);
    for (CellIterator cit = cells.begin();
        cit != cells.end();
        ++cit)
    {

      const long row_index2 = lattice_temperature.get_unknown_index(*cit);
      if (row_index2 < 0)  //here is a Dirichlet boundary condition
        continue;
      const std::size_t row_index = std::size_t(row_index2);

      PointType centroid_cell = viennagrid::centroid(*cit);

      const double T_center = lattice_temperature.get_value(*cit);

      const double kappa_center = diffusivity(*cit, T_center);

      b[row_index]   = 0;
      A(row_index, row_index) = 0;

      //
      //   K * laplace T
      //
      FacetOnCellContainer facets(*cit);
      for (FacetOnCellIterator focit = facets.begin();
          focit != facets.end();
          ++focit)
      {
        CellType const *other_cell_ptr = util::get_other_cell_of_facet(mesh, *focit, *cit);

        if (!other_cell_ptr) continue;

        const long col_index  = lattice_temperature.get_unknown_index(*other_cell_ptr);

        PointType centroid_other_cell = viennagrid::centroid(*other_cell_ptr);
        PointType cell_connection = centroid_other_cell - centroid_cell;
        PointType cell_connection_normalized = cell_connection / viennagrid::norm(cell_connection);
        PointType facet_unit_normal = viennashe::util::outer_cell_normal_at_facet(*cit, *focit);

        const double connection_len = viennagrid::norm_2(centroid_cell - centroid_other_cell);
        const double weighted_interface_area = viennagrid::volume(*focit) * viennagrid::inner_prod(facet_unit_normal, cell_connection_normalized);

        const double T_outer = lattice_temperature.get_value(*other_cell_ptr);

        // off-diagonal contribution
        if (col_index >= 0)
        {
          PointType facet_center = viennagrid::centroid(*focit);  //TODO: Use intersection of facet plane with connection
          const double connection_in_cell = viennagrid::norm(facet_center - centroid_cell);
          const double connection_in_other_cell = viennagrid::norm(facet_center - centroid_other_cell);
          const double kappa_mean = (connection_in_cell + connection_in_other_cell) /
                                       (connection_in_cell/kappa_center + connection_in_other_cell/diffusivity(*other_cell_ptr, T_outer));

          A(row_index, std::size_t(col_index))  = kappa_mean * weighted_interface_area / connection_len;
          A(row_index, row_index) -= kappa_mean * weighted_interface_area / connection_len;

          b[row_index] -= kappa_mean * weighted_interface_area * (T_outer - T_center) / connection_len;
        }
        else if (lattice_temperature.get_boundary_type(*other_cell_ptr) == BOUNDARY_DIRICHLET)
        {
          A(row_index, row_index) -= kappa_center * weighted_interface_area / connection_len;

          b[row_index] -= kappa_center * weighted_interface_area / connection_len * (lattice_temperature.get_boundary_value(*other_cell_ptr) - T_center);
        }

      } //for facets

      const double volume = viennagrid::volume(*cit);

      // residual contribution
      b[row_index] += volume * quan_power_density(*cit); // / kappa_center

    } //for vertices

  } //assemble_poisson

  /**
   * @brief Assemble a spatial quantity, where the equation is deduced from 'quan'.
   * @param device The device ...
   * @param quantities A list of required and known quantities
   * @param conf The simulator configuration
   * @param quan The unkown quantity
   * @param A The system matrix
   * @param b The right hand side
   *
   * @throw Throws an assembly_not_implemented_exception if no assembly routine for quan can be found.
   *
   */
  template <typename DeviceType,
            typename TimeStepQuantitiesT,
            typename VertexT,
            typename MatrixType,
            typename VectorType>
  void assemble( DeviceType const & device,
                 TimeStepQuantitiesT & quantities,
                 viennashe::config const & conf,
                 viennashe::unknown_quantity<VertexT> const & quan,
                 MatrixType & A,
                 VectorType & b)
  {
    //
    // Improve this by using ViennaMath to dispatch
    //

    if (quan.get_name() == viennashe::quantity::potential())
      assemble_poisson(device, quantities, conf, A, b);
    else if (quan.get_name() == viennashe::quantity::electron_density())
      assemble_dd(device, quantities, conf, ELECTRON_TYPE_ID, A, b);
    else if (quan.get_name() == viennashe::quantity::hole_density())
      assemble_dd(device, quantities, conf,     HOLE_TYPE_ID, A, b);
    else if (quan.get_name() == viennashe::quantity::density_gradient_electron_correction())
      assemble_density_gradient(device, quantities, conf, ELECTRON_TYPE_ID, A, b);
    else if (quan.get_name() == viennashe::quantity::density_gradient_hole_correction())
      assemble_density_gradient(device, quantities, conf,     HOLE_TYPE_ID, A, b);
    else if (quan.get_name() == viennashe::quantity::lattice_temperature())
      assemble_heat(device, quantities, conf, A, b);
    else
      throw assembly_not_implemented_exception(quan.get_name());
  }

} //namespace viennashe

#endif
