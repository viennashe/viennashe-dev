#ifndef VIENNASHE_SHE_SIMULATOR_HPP
#define VIENNASHE_SHE_SIMULATOR_HPP

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

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/device.hpp"
#include "viennashe/simulator_setup.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/assemble.hpp"

#include "viennashe/mapping.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/log_keys.h"

#include "viennashe/math/linalg_util.hpp"

#include "viennashe/she/assemble_all.hpp"
#include "viennashe/she/expansion_order.hpp"
#include "viennashe/she/eliminate.hpp"
#include "viennashe/she/scattering/all.hpp"
#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/she/setup_energies.hpp"
#include "viennashe/she/linear_solver.hpp"
#include "viennashe/she/setup_energies.hpp"
#include "viennashe/she/timestep_quantities.hpp"
#include "viennashe/she/df_wrappers.hpp"
#include "viennashe/she/boundary_conditions.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/she/transfer_to_new_h_space.hpp"
#include "viennashe/she/log_keys.h"

#include "viennashe/util/timer.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"
#include "viennashe/util/smooth_doping.hpp"

// HDE
#include "viennashe/phonon/joule_heating.hpp"

/** @file viennashe/simulator.hpp
    @brief Implements the SHE simulator classes (both self-consistent and non-self-consistent).
*/

namespace viennashe
{

  namespace detail
  {

    /**
     * @brief Sets the unknown masks per material.
     *
     * @param device The device for which to set the unkown masks
     * @param quan The unkown quantity for which to set the mask
     * @param material_check A functor, which takes a ViennaGrid cell and
     *                       returns true if the mask has to be set on the given cell (does a material check)
     */
    template<typename DeviceT, typename QuantityT>
    void set_unknown_for_material(DeviceT const & device, QuantityT & quan, materials::checker const & material_check)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (!material_check(device.get_material(*cit)))
          continue;

        quan.set_unknown_mask(*cit, true);
      } // for cells
    } // set_unknown_for_material

    template<typename DeviceT, typename QuantityT>
    void set_she_unknown(DeviceT const & device, QuantityT & quan)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        quan.set_unknown_mask(*cit, false); //initialize

        if (viennashe::materials::is_semiconductor(device.get_material(*cit))) // Semiconductor: good
          quan.set_unknown_mask(*cit, true);
        else if (viennashe::materials::is_conductor(device.get_material(*cit))) // Conductor: Set unknown if connected to semiconductor
        {
          viennagrid_element_id *semi_cell = viennashe::util::get_connected_semiconductor_cell(device, *cit);

          if (semi_cell)
            quan.set_unknown_mask(*cit, true);
        }
      } // for cells
    } // set_unknown_for_material


    /**
     * @brief Sets the boundary condition (type and value) per vertex depending on the material
     *
     * @param device The device for which to set boundary conditions
     * @param quan The unkown quantity for which to set boundary conditions
     * @param material_check A functor, which is given a ViennaGrid cell and
     *                       returns true if there are boundary conditions to be set on the vertices of that cell
     * @param boundary_value_accessor An accessor, which returns the boundary value given a ViennaGrid vertex
     * @param bnd_id The id of the boundary type
     */
    template<typename DeviceT, typename QuantityT, typename BoundaryValueAccessorT>
    void set_boundary_for_material(DeviceT const & device,
                                   QuantityT & quan,
                                   materials::checker const & material_check,
                                   BoundaryValueAccessorT boundary_value_accessor,
                                   boundary_type_id bnd_id)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (!material_check(device.get_material(*cit)))
          continue;

        typename BoundaryValueAccessorT::value_type bnd_value = boundary_value_accessor(*cit);
        quan.set_boundary_type(*cit, bnd_id);
        quan.set_boundary_value(*cit, bnd_value);
        quan.set_value(*cit, bnd_value);
      } // for cells
    } // set_boundary_for_material


    /**
     * @brief Sets the boundary condition type per material on vertices. Does not set the boundary value
     *
     * @param device The device
     * @param quan The unkown quantity for which to set the boundary condition type
     * @param material_check A functor, which is given a ViennaGrid cell and
     *                       returns true if there are boundary conditions to be set on the vertices of that cell
     * @param bnd_id The id of the boundary type
     */
    template<typename DeviceT>
    void set_boundary_for_material(DeviceT const & device,
                                   viennashe::she::unknown_she_quantity<double> & quan,
                                   materials::checker const & material_check,
                                   boundary_type_id bnd_id)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        if (!material_check(device.get_material(*cit)))
          continue;

        quan.set_boundary_type(*cit, bnd_id);
      }
    }

    /**
     * @brief Sets the initial guess per quantity and device
     * @param device The device
     * @param quan The unkown quantity for which the initial guess has to be set
     * @param initial_guess_accessor An accessor, which returns a value given a ViennaGrid vertex
     */
    template<typename DeviceT, typename QuantityT, typename BoundaryValueAccessorT>
    void set_initial_guess(DeviceT const & device,
                           QuantityT & quan,
                           BoundaryValueAccessorT const & initial_guess_accessor)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        quan.set_value(*cit, initial_guess_accessor(*cit));
      }
    }

  } // namespace detail


  // TODO: Think about whether this is a useful way of dealing with temperature. Doesn't seem right.
  template<typename DeviceT, typename UnknownQuantityListT>
  void transfer_quantity_to_device(DeviceT & device, UnknownQuantityListT const & quantities, viennashe::config const & conf)
  {
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    //
    // Transfer lattice temperature
    if (conf.with_hde())
    {
      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        const double TL = quantities.get_unknown_quantity(viennashe::quantity::lattice_temperature()).get_value(*cit);
        device.set_lattice_temperature(TL, *cit);
      }
    }

    //
    // Insert transfer for other quantities here
    //

  } // transfer_quantity_to_device

  /**
   * @brief Returns the relative update norm (L_2 ; realtive to the current value). Except for the potential.
   * @param spatial_quan The spatial quantity
   * @param x The solution vector
   * @return The relative update norm (L_2)
   */
  template<typename UnknownQuantityT, typename VectorType>
  double get_relative_update_norm(UnknownQuantityT const & spatial_quan, VectorType const & x)
  {

    double norm = 0;
    for ( std::size_t j = 0; j < x.size(); ++j)
    {
      if (spatial_quan.get_name() == viennashe::quantity::potential())
      {
        norm += x.at(j) * x.at(j);
      }
      else
      {
        const double current_val = spatial_quan.values().at(j);
        const double zw          = (std::abs(current_val) > 0) ? (x.at(j) / current_val) : ( x.at(j) );
        norm += zw * zw ;
      }
    }
    return sqrt(norm);
  } // get_relative_update_norm

  /**
   * @brief Returns the relative update norm for SHE based carrier concentrations (L_2 norm).
   * @param device The device
   * @param quan The unknown SHE quantity
   * @param spatial_quan The spatial quantity to update (e.g. n or p)
   * @param conf The simulator configuration
   * @return An relative update norm for SHE based spatial quantities
   */
  template<typename DeviceT>
  double get_relative_update_norm(DeviceT const & device,
                       viennashe::she::unknown_she_quantity<double> const & quan,
                       viennashe::unknown_quantity<double> const & spatial_quan,
                       viennashe::config const & conf)
  {
    typedef viennashe::she::unknown_she_quantity<double>   SHEQuantityType;

    viennashe::she::carrier_density_wrapper<SHEQuantityType> density_wrapper(conf, quan);

    double norm = 0;

    // Set bandedge shift on vertices and edges:
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    std::vector<double> tmp(3);
    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      tmp = density_wrapper(*cit);
      const double new_value = tmp[0];
      const double old_value = spatial_quan.get_value(*cit);
      double zw  = 0;
      if (std::fabs(old_value) > 0)
        zw = (new_value - old_value) / old_value;
      else
        zw = new_value;

      norm += zw * zw;
    }

    return std::sqrt(norm);
  } // get_relative_update_norm



  /**
   * @brief Potential update calculator for Newton.
   * @param device               The device
   * @param unknown_quantities   The unkown quantity to be processed
   * @param x                    The update vector for all quantities
   */
  template<typename DeviceT, typename TimeStepQuantitiesT, typename VectorT>
  double get_potential_update_norm_newton(DeviceT const & device,
                                          TimeStepQuantitiesT const & unknown_quantities,
                                          VectorT const & x)
  {
    typedef typename TimeStepQuantitiesT::UnknownQuantityType   UnknownQuantityType;

    double potential_update_norm = 0;

    for (std::size_t i=0; i<unknown_quantities.unknown_quantities().size(); ++i)
    {
      UnknownQuantityType const & current_quan = unknown_quantities.unknown_quantities()[i];

      if (current_quan.get_name() != viennashe::quantity::potential())
        continue;

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        long idx = current_quan.get_unknown_index(*cit);
        if (idx >= 0)
          potential_update_norm += x[std::size_t(idx)] * x[std::size_t(idx)];
      }
    }

    return std::sqrt(potential_update_norm);
  } // get_potential_update_norm_newton

  /**
   * @brief Total update calculator for Newton. Needs to be called *before* update_quantities!
   * @param device The device
   * @param unknown_quantities_ The unkown quantity to be processed
   * @param x The update to the quantity
   */
  template<typename DeviceT, typename UnknownQuantityT, typename VectorT>
  double get_total_update_norm_newton(DeviceT const & device,
                                      UnknownQuantityT const & unknown_quantities_,
                                      VectorT const & x)
  {
    (void)device; (void)unknown_quantities_;
    double total_update_norm = 0;

    for ( std::size_t j = 0; j < x.size(); ++j)
        total_update_norm += x.at(j) * x.at(j);

    return total_update_norm;
  } // get_total_update_norm_newton

  /**
   * @brief Updater for Newton
   * @param device             The device
   * @param conf               The simulator configuration (cf. config.hpp)
   * @param unknown_quantities The unkown quantity to be updated
   * @param x                  The vector containing the updates to the quantity
   */
  template<typename DeviceT, typename TimeStepQuantitiesT, typename VectorType>
  void update_quantities(DeviceT const & device,
                         viennashe::config const & conf,
                         TimeStepQuantitiesT & unknown_quantities,
                         VectorType const & x)
  {
    typedef typename TimeStepQuantitiesT::UnknownQuantityType   UnknownQuantityType;

    //
    // Step 1: Prepare update and current value arrays (split up into the various quantities)
    //
    for (std::size_t i=0; i<unknown_quantities.unknown_quantities().size(); ++i)
    {
      UnknownQuantityType & current_quan = unknown_quantities.unknown_quantities()[i];

      double norm_inf_update  = 0;
      double norm_inf_current = 0;

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        long idx = current_quan.get_unknown_index(*cit);
        if (idx >= 0)
        {
          norm_inf_current = std::max(norm_inf_current, current_quan.get_value(*cit));
          norm_inf_update  = std::max(norm_inf_update,  x[std::size_t(idx)]);
        }
      }

      double norm_rel_update = norm_inf_current > 0 ? norm_inf_update / norm_inf_current : 1.0;

      if (norm_inf_update == 0.0
          //|| norm_rel_update < 1e-14 //no significant update at all, so move on to next quantity
         )
        continue;

      double delta = 0.1;

      std::size_t damping_iter_max = 10;
      double alpha = 1.0;

      for (std::size_t damping_iter = 0; damping_iter < damping_iter_max; ++damping_iter)
      {
        if (damping_iter == 0) //check whether a full Newton step is okay
          alpha = 1.0;
        else if (damping_iter == 1) //use logarithmic damping, cf. PhD thesis by Simlinger or Fischer
          alpha = (1.0 + delta * std::log(norm_rel_update)) / (1.0 + delta * (norm_rel_update - 1.0));
        else
          alpha = conf.nonlinear_solver().damping() * alpha;

        double update_norm = 0.0;
        double current_norm = 0.0;

        bool do_reject = false;

        // Check to see whether new quantity makes sense:
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          long index = current_quan.get_unknown_index(*cit);
          if (index >= 0)
          {
            double  update_value = alpha * x[std::size_t(index)];
            double current_value = current_quan.get_value(*cit);

            update_norm  +=  update_value *  update_value;
            current_norm += current_value * current_value;

            if (   std::abs(current_value) != 0.0
                && std::abs(current_value) < 0.5 * std::abs(alpha * update_value)
               )
              do_reject = true;

            if (current_quan.get_logarithmic_damping() && current_value + update_value <= 0)
              do_reject = true;
          }
        }

        if (update_norm > current_norm)
          do_reject = true;

        // Don't accept updated variables if update is too large
        if ( (damping_iter < damping_iter_max - 1) //last iterate is always accepted
            && do_reject)
        {
          continue;
        }

        // write new values:
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          long index = current_quan.get_unknown_index(*cit);
          if (index >= 0)
          {
            double  update_value = alpha * x[std::size_t(index)];
            double current_value = current_quan.get_value(*cit);

            current_quan.set_value(*cit, current_value + update_value);
          }
        }
        break;

      } //for damping_iter


    } //for quantities
  }

  /** @brief Updater for SHE quantities when using Gummel iteration
   *
   * @param device The device
   * @param quan The SHE unkown quantity to be updated
   * @param spatial_quan The spatial, eg. n or p or potential, unkown quantity
   * @param conf The simulator configuration (cf. config.hpp)
   * @param x The update to the SHE quantity
   * @param force_no_damping If true no damping will be applied, that is d = 1.0
   */
  template<typename DeviceT, typename VectorType>
  void update_quantity(DeviceT const & device,
                       viennashe::she::unknown_she_quantity<double> & quan,
                       viennashe::unknown_quantity<double> & spatial_quan,
                       viennashe::config const & conf,
                       VectorType const & x, bool force_no_damping = false)
  {
    typedef viennashe::she::unknown_she_quantity<double>   SHEQuantityType;

    typedef viennashe::she::detail::carrier_density_wrapper_by_reference<SHEQuantityType> density_wrapper_type;

    density_wrapper_type density_wrapper(conf, quan);

    // Set bandedge shift on vertices and edges:
    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
      {
        long index = quan.get_unknown_index(*cit, index_H);
        if (index >= 0)
          quan.set_values(*cit, index_H, &(x[static_cast<std::size_t>(index)]));
      }

      // write density:
      if (force_no_damping)
      {
        spatial_quan.set_value(*cit, density_wrapper(*cit));
      }
      else
      {
        const double alpha         = conf.nonlinear_solver().damping();
              double new_value     = density_wrapper(*cit);
        const double current_value = spatial_quan.get_value(*cit);
        if (spatial_quan.get_logarithmic_damping())
        {
          if (new_value < current_value * 1e-10) // prevent problems with fractional powers of negative numbers by limiting the update
            new_value = current_value * 1e-10;

          const double a1 = std::pow(current_value, 1.0 - alpha);
          const double a2 = std::pow(new_value, alpha);

          spatial_quan.set_value(*cit, a1 * a2);
        }
        else
        {
          spatial_quan.set_value(*cit, current_value * (1.0 - alpha) + new_value * alpha);
        }
      }
    } // for vertices

    viennagrid_element_id *facets_begin, *facets_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim - 1, &facets_begin, &facets_end));
    for (viennagrid_element_id *fit  = facets_begin;
                                fit != facets_end;
                              ++fit)
    {
      for (std::size_t index_H = 0; index_H < quan.get_value_H_size(); ++index_H)
      {
        long index = quan.get_unknown_index(*fit, index_H);
        if (index >= 0)
          quan.set_values(*fit, index_H, &(x[static_cast<std::size_t>(index)]));
      }
    }

  } // update_quantity

  /** @brief Updater for spatial unknowns (i.e. defined in x-space) when using Gummel iteration
   *
   * @param device The device
   * @param unknown_quantity The unkown quantity to be updated
   * @param alpha The damping, bounded by (0;1]
   * @param x The update to the quantity
   */
  template<typename DeviceT, typename VectorType>
  void update_quantity(DeviceT const & device,
                       viennashe::unknown_quantity<double> & unknown_quantity,
                       double alpha,
                       VectorType const & x)
  {
    double update_norm = 0;

    viennagrid_dimension cell_dim;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

    viennagrid_element_id *cells_begin, *cells_end;
    VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      double current_value = unknown_quantity.get_value(*cit);
      double new_value     = current_value;

      if (unknown_quantity.get_boundary_type(*cit) == BOUNDARY_DIRICHLET) //boundary value
        new_value = unknown_quantity.get_boundary_value(*cit);
      else if (unknown_quantity.get_unknown_mask(*cit)) // interior value
        new_value = current_value + x[std::size_t(unknown_quantity.get_unknown_index(*cit))];

      update_norm += std::abs(new_value - current_value);

      if (unknown_quantity.get_logarithmic_damping())
      {
        if (new_value < current_value * 1e-10) // prevent problems with fractional powers of negative numbers by limiting the update
        {
          new_value = current_value * 1e-10;
        }

        unknown_quantity.set_value(*cit, std::pow(current_value, 1.0 - alpha) * std::pow(new_value, alpha));
      }
      else
      {
        unknown_quantity.set_value(*cit, current_value * (1.0 - alpha) + new_value * alpha);
      }
    }

    //log::info<log_simulator>() << "Norm of update (undamped): " << update_norm << std::endl;
  } // update_quantity

  /** @brief Updater for spatial unknowns (i.e. defined in x-space) when using Gummel iteration
   *
   * @param device The device
   * @param unknown_quantity The unkown quantity to be updated
   * @param conf The simulator configuration (cf. config.hpp)
   * @param x The update to the quantity
   * @param force_no_damping If true no damping will be applied, that is d = 1.0
   */
  template<typename DeviceT, typename VectorType>
  void update_quantity(DeviceT const & device,
                       viennashe::unknown_quantity<double> & unknown_quantity,
                       viennashe::config const & conf,
                       VectorType const & x, bool force_no_damping = false)
  {
    if (force_no_damping) update_quantity(device, unknown_quantity, 1.0, x);
    else                  update_quantity(device, unknown_quantity, conf.nonlinear_solver().damping(), x);
  }

  /** @brief  Class for self-consistent SHE simulations.
   *
   * @tparam DeviceType      Type of the device the simulator is operating on
   */
  template <typename DeviceType>
  class simulator
  {
    private:
      typedef simulator<DeviceType>    self_type;

      typedef viennashe::math::sparse_matrix<double>    MatrixType;
      typedef std::vector<double>                       VectorType;

      typedef typename DeviceType::mesh_type           MeshType;

    public:

      typedef DeviceType device_type;

      typedef viennashe::she::timestep_quantities<DeviceType>            SHETimeStepQuantitiesT;

      typedef typename SHETimeStepQuantitiesT::she_df_type                       she_df_type;
      typedef typename SHETimeStepQuantitiesT::edf_type                             edf_type;
      typedef typename SHETimeStepQuantitiesT::generalized_edf_type     generalized_edf_type;

      typedef typename SHETimeStepQuantitiesT::UnknownSHEQuantityType      she_quantity_type;

      typedef unknown_quantity<double>       UnknownQuantityType;
      typedef UnknownQuantityType            unknown_quantity_type;

      typedef const_quantity<double>         ResultQuantityType;

      typedef ResultQuantityType          potential_type;
      typedef ResultQuantityType   electron_density_type;
      typedef ResultQuantityType       hole_density_type;


      /** @brief Constructs the self-consistent simulator object
       *
       * @param device  The device
       * @param conf    A SHE configuation object
       */
      simulator(DeviceType & device, viennashe::config const & conf = viennashe::config()) : p_device_(&device), config_(conf)
      {
        quantities_history_.push_back(SHETimeStepQuantitiesT());

        // Step 1: Doping must be available on cells
        //setup_doping_on_vertices(device);

        if(conf.setup_insulator_distances())
        {
          throw std::runtime_error("TODO: setup_insulator_distances()");
          //setup_insulator_distances(device);
        }

        // ensure doping in vicinity of contact is constant
        detail::smooth_doping_at_contacts(device);

        // push quantities
        // (note that by default all values are 'known' default values, so one has to specify the unknown regions later):
        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

        viennagrid_int  cell_num;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(device.mesh(), cell_dim, &cell_num));
        viennagrid_int facet_num;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(device.mesh(), cell_dim - 1, &facet_num));

        // electrostatic potential:
        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::potential(), EQUATION_POISSON_DD, cell_num));
        detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_NO_CONDUCTOR_ID));
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), boundary_potential_accessor<DeviceType>(device), BOUNDARY_DIRICHLET);
        detail::set_initial_guess(device, quantities().unknown_quantities().back(), built_in_potential_accessor<DeviceType>(device));

        // electron and hole densities:
        contact_carrier_density_accessor<DeviceType> contact_carrier_density_holes(device, HOLE_TYPE_ID, true);
        contact_carrier_density_accessor<DeviceType> contact_carrier_density_electrons(device, ELECTRON_TYPE_ID, true);
        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::electron_density(), EQUATION_CONTINUITY, cell_num, 1.0));
        if (conf.with_electrons())
        {
          if (conf.get_electron_equation() == EQUATION_CONTINUITY) //enable standard DD instead of SHE
            detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_SEMICONDUCTOR_ID));
        }
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), contact_carrier_density_electrons, BOUNDARY_DIRICHLET);

        detail::set_initial_guess(device, quantities().unknown_quantities().back(), contact_carrier_density_electrons);

        quantities().unknown_quantities().back().set_logarithmic_damping(true);

        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::hole_density(), EQUATION_CONTINUITY, cell_num, 1.0));
        if (conf.with_holes())
        {
          if (conf.get_hole_equation() == EQUATION_CONTINUITY) // enable standard DD instead of SHE
            detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_SEMICONDUCTOR_ID));
        }
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), contact_carrier_density_holes, BOUNDARY_DIRICHLET);


        detail::set_initial_guess(device, quantities().unknown_quantities().back(), contact_carrier_density_holes);

        quantities().unknown_quantities().back().set_logarithmic_damping(true);

        // Density gradient: Correction potential for ELECTRONS
        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::density_gradient_electron_correction(), EQUATION_DENSITY_GRADIENT, cell_num, 0.0));
        if (conf.quantum_correction()) detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_SEMICONDUCTOR_ID));
        // BND Dirichlet with gamma=0 at conductors
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), constant_accessor<double>(0), BOUNDARY_DIRICHLET);
        // BND for insulators
        if (conf.density_gradient(viennashe::ELECTRON_TYPE_ID).boundary_type() == viennashe::BOUNDARY_DIRICHLET)
        {
          if (conf.quantum_correction())
            viennashe::log::info<viennashe::log_simulator>() << "* simulator(): Density Gradient (electrons) ... Insulator boundary condition is DIRICHLET " << std::endl;
          const double dg_bnd_value = conf.density_gradient(viennashe::ELECTRON_TYPE_ID).dirichlet_boundary_value();
          detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_INSULATOR_ID), constant_accessor<double>(dg_bnd_value), BOUNDARY_DIRICHLET);
        }
        else
        {
          if (conf.quantum_correction())
            viennashe::log::info<viennashe::log_simulator>() << "* simulator(): Density Gradient (electrons) ... Insulator boundary condition is ROBIN " << std::endl;
          robin_boundary_coefficients<double> robin_coeffs_n = conf.density_gradient(viennashe::ELECTRON_TYPE_ID).robin_coeffs();
          detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_INSULATOR_ID), constant_accessor<robin_boundary_coefficients<double> >(robin_coeffs_n), BOUNDARY_ROBIN);
        }

        // Density gradient: Correction potential for HOLES
        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::density_gradient_hole_correction(),     EQUATION_DENSITY_GRADIENT, cell_num, 0.0));
        if (conf.quantum_correction()) detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_SEMICONDUCTOR_ID));
        // BND Dirichlet with gamma=0 at conductors
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), constant_accessor<double>(0), BOUNDARY_DIRICHLET);
        // BND for insulators
        if (conf.density_gradient(viennashe::HOLE_TYPE_ID).boundary_type() == viennashe::BOUNDARY_DIRICHLET)
        {
          if (conf.quantum_correction())
            viennashe::log::info<viennashe::log_simulator>() << "* simulator(): Density Gradient (holes) ... Insulator boundary condition is DIRICHLET " << std::endl;
          const double dg_bnd_value = conf.density_gradient(viennashe::HOLE_TYPE_ID).dirichlet_boundary_value();
          detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_INSULATOR_ID), constant_accessor<double>(dg_bnd_value), BOUNDARY_DIRICHLET);
        }
        else
        {
          if (conf.quantum_correction())
            viennashe::log::info<viennashe::log_simulator>() << "* simulator(): Density Gradient (holes) ... Insulator boundary condition is ROBIN " << std::endl;
          robin_boundary_coefficients<double> robin_coeffs_p = conf.density_gradient(viennashe::HOLE_TYPE_ID).robin_coeffs();
          detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_INSULATOR_ID), constant_accessor<robin_boundary_coefficients<double> >(robin_coeffs_p), BOUNDARY_ROBIN);
        }

        // Temperature
        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::lattice_temperature(), EQUATION_POISSON_HEAT,  cell_num, 300.0));
        if (conf.with_hde())
          detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_NO_CONDUCTOR_ID));
        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID),
                                          viennashe::lattice_temperature_accessor<DeviceType>(device), BOUNDARY_DIRICHLET);
        // initial guess already set above

        //
        // Traps
//        quantities().unknown_quantities().push_back(UnknownQuantityType(viennashe::quantity::trap_occupancy(), viennashe::EQUATION_SRH_TRAPPING,  cell_num, 0.0));
//        if (conf.with_traps())
//        {
//          // At the moment we only support bulk SRH traps
//          detail::set_unknown_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_SEMICONDUCTOR_ID));
//        }
//        detail::set_boundary_for_material(device, quantities().unknown_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID),
//                                          constant_accessor<double>(0.0), BOUNDARY_DIRICHLET);

        quantities().setup_trap_unkown_indices(this->device());


        // SHE: electron distribution function
        quantities().unknown_she_quantities().push_back(she_quantity_type(p_device_->mesh(), cell_dim, cell_dim - 1, viennashe::quantity::electron_distribution_function(), ELECTRON_TYPE_ID, EQUATION_SHE));
        quantities().unknown_she_quantities().back().resize(cell_num, facet_num);
        if (conf.with_electrons())
        {
          if (conf.get_electron_equation() == EQUATION_SHE) //enable standard DD instead of SHE
            detail::set_she_unknown(device, quantities().unknown_she_quantities().back());
        }
        detail::set_boundary_for_material(device, quantities().unknown_she_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), BOUNDARY_DIRICHLET);
        quantities().unknown_quantities().back().set_logarithmic_damping(true);

        // SHE: hole distribution function
        quantities().unknown_she_quantities().push_back(she_quantity_type(p_device_->mesh(), cell_dim, cell_dim - 1, viennashe::quantity::hole_distribution_function(),         HOLE_TYPE_ID, EQUATION_SHE));
        quantities().unknown_she_quantities().back().resize(cell_num, facet_num);
        if (conf.with_holes())
        {
          if (conf.get_hole_equation() == EQUATION_SHE) //enable standard DD instead of SHE
            detail::set_she_unknown(device, quantities().unknown_she_quantities().back());
        }
        detail::set_boundary_for_material(device, quantities().unknown_she_quantities().back(), materials::checker(MATERIAL_CONDUCTOR_ID), BOUNDARY_DIRICHLET);
      }

      /**
       * @brief Transfers the inital guess for the given quantity
       *
       * @param quan_name The name of the quantity
       * @param quan_acc An accessor, which is given a ViennaGrid vertex and returns a value
       */
      template <typename QuantityAccessorT>
      void set_initial_guess(std::string quan_name, QuantityAccessorT const & quan_acc)
      {
        transfer_provided_quantities(quan_acc, quantities().get_unknown_quantity(quan_name));
      }


      /** @brief Launches the solver. Uses the built-in potential as initial guess for the potential and
       *         the doping concentration as the initial guess for carriers.
       */
      void run()
      {
        const double use_newton = (config().nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

        detail::set_boundary_for_material(device(), quantities().unknown_she_quantities()[0], materials::checker(MATERIAL_CONDUCTOR_ID), BOUNDARY_DIRICHLET);
        detail::set_boundary_for_material(device(), quantities().unknown_she_quantities()[1], materials::checker(MATERIAL_CONDUCTOR_ID), BOUNDARY_DIRICHLET);

        //
        // Run nonlinear iterations
        //

        double initial_residual_norm = 0;

        SHETimeStepQuantitiesT transferred_quantities;

        bool update_no_damping  = false;
        bool break_after_update = false;

        for (std::size_t nonlinear_iter = 1; nonlinear_iter <= this->config().nonlinear_solver().max_iters(); ++nonlinear_iter)
        {
          viennashe::util::timer stopwatch;
          stopwatch.start();

          // If this is the last iteration set update_no_damping
          if (nonlinear_iter == this->config().nonlinear_solver().max_iters())
            update_no_damping = true;

          //
          // write kinetic energy to device:
          //
          //log::debug<viennashe::she::log_she_solver>() << "* simulator(): Writing center of band gap to device..." << std::endl;
          viennashe::she::setup_energies(this->device(), this->quantities(), this->config());

          //
          // distribute SHE expansion orders over device
          //
          //log::debug<viennashe::she::log_she_solver>() << "* simulator(): Writing expansion orders to device (" << config().max_expansion_order() << ")... " << std::endl;
          viennashe::she::distribute_expansion_orders(this->device(), this->quantities(), this->config());

          //
          // write boundary conditions:
          //
          //log::debug<viennashe::she::log_she_solver>() << "* simulator(): Writing boundary conditions to device..." << std::endl;
          viennashe::she::write_boundary_conditions(this->device(), this->quantities(), this->config());

          //
          // Mapping of unknowns:
          //
          viennashe::map_info_type map_info = viennashe::create_mapping(this->device(), this->quantities(), this->config());

          transferred_quantities = quantities(); //transfer kinetic energy and other information related to the (x,H)-space
          if ( quantities_history_.size() > 1 )
          {
            transferred_quantities.get_unknown_quantity(viennashe::quantity::potential()) = quantities_history_.at(quantities_history_.size()-2).get_unknown_quantity(viennashe::quantity::potential());
            viennashe::she::transfer_to_new_h_space(device(), quantities_history_.at(quantities_history_.size()-2), transferred_quantities, this->config());

            // push every nonlinear iteration state:
            //quantities_history_.at(quantities_history_.size()-2) = transferred_quantities;
          }

          //
          // Transfer device based quantities back to the device
          //
          viennashe::transfer_quantity_to_device(this->device(), this->quantities(), this->config());


          if (nonlinear_iter == 1)
          {
            for (std::size_t i = 0; i < this->quantities().unknown_quantities().size(); ++i)
            {
              log::info() << "* run(): Number of unknowns for quantity '" << this->quantities().unknown_quantities()[i].get_name() << "': "
                          << map_info[this->quantities().unknown_quantities()[i].get_name()].first << std::endl;
            }
          }

          // Log summary output
          if (log_simulator::enabled)
          {
            if (nonlinear_iter == 1)
            {
              //
              // Write convergence table header:
              //
              if (use_newton)
                log::info<log_simulator>() << "* run(): -- Newton scheme --" << std::endl;
              else
                log::info<log_simulator>() << "* run(): -- Gummel scheme --" << std::endl;

              log::info<log_simulator>() << std::setw(5) <<  " Iter" << " | ";

              if (map_info[viennashe::quantity::electron_distribution_function()].first > 0)
                log::info<log_simulator>() << std::setw(9) << "#f_even^e" << " | " << std::setw(8) << "#f_odd^e" << " | ";
              if (map_info[viennashe::quantity::hole_distribution_function()].first > 0)
                log::info<log_simulator>() << std::setw(8) << "#f_even^h" << " | " << std::setw(8) << "#f_odd^h" << " | ";

              log::info<log_simulator>() << std::setw(9) << "ResidNorm"      << " | "
                                         << std::setw(9) << "|d_pot|_2"      << " | "
                                         << std::setw(9) << "TotUpdate"      << " | "
                                         << std::setw(8) << "Time (s)"     << std::endl;
            }

            log::info<log_simulator>() << std::setw(5) << nonlinear_iter << " | ";

            if (map_info[viennashe::quantity::electron_distribution_function()].first > 0)
              log::info<log_simulator>() << std::setw(9) << viennashe::util::format_number_to_width(map_info[viennashe::quantity::electron_distribution_function()].first, 7) << " | "
                                         << std::setw(8) << viennashe::util::format_number_to_width(map_info[viennashe::quantity::electron_distribution_function()].second, 7) << " | ";
            if (map_info[viennashe::quantity::hole_distribution_function()].first > 0)
              log::info<log_simulator>() << std::setw(9) << viennashe::util::format_number_to_width(map_info[viennashe::quantity::hole_distribution_function()].first, 7) << " | "
                                         << std::setw(8) << viennashe::util::format_number_to_width(map_info[viennashe::quantity::hole_distribution_function()].second, 7) << " | ";

            std::cout << std::flush ;
          }

          //
          // Nonlinear methods (Newton, Gummel)
          //

          double potential_norm_increment = 0;
          double current_residual_norm    = 0;
          double total_update_norm        = 0;

          if (use_newton)
          {
            std::size_t total_number_of_unknowns = 0;
            for (viennashe::map_info_type::const_iterator it = map_info.begin(); it != map_info.end(); ++it)
              total_number_of_unknowns += (it->second).first + (it->second).second;   //sum of unknowns on vertices and edges

            MatrixType A(total_number_of_unknowns, total_number_of_unknowns);
            VectorType b(total_number_of_unknowns);

            // assemble spatial quantities:
            for (std::size_t i = 0; i < this->quantities().unknown_quantities().size(); ++i)
              viennashe::assemble(device(), this->quantities(), this->config(), this->quantities().unknown_quantities()[i], A, b);

            // assemble SHE quantities:
            for (std::size_t i = 0; i < this->quantities().unknown_she_quantities().size(); ++i)
              viennashe::she::assemble(this->device(), transferred_quantities, this->quantities(), this->config(), this->quantities().unknown_she_quantities()[i], A, b,
                                        (quantities_history_.size() > 1), nonlinear_iter > 1);

            //VectorType x = viennashe::she::solve(A, b, map_info);  //TODO: use this!
            VectorType x = this->solve(A, b, true);

            current_residual_norm += viennashe::math::norm_2(b);

            total_update_norm        = viennashe::get_total_update_norm_newton(this->device(), this->quantities(), x);
            potential_norm_increment = viennashe::get_potential_update_norm_newton(this->device(), this->quantities(), x);

            viennashe::update_quantities(this->device(), this->config(), this->quantities(), x);
          }
          else // Gummel iteration
          {
            // spatial quantities:
            for (std::size_t i=0; i < this->quantities().unknown_quantities().size(); ++i)
            {
              std::size_t number_of_unknowns = map_info[quantities().unknown_quantities()[i].get_name()].first;
              if (number_of_unknowns == 0)
                continue;

              // System for this quantity only:
              MatrixType A(number_of_unknowns, number_of_unknowns);
              VectorType b(number_of_unknowns);

              viennashe::assemble(this->device(), this->quantities(), this->config(), this->quantities().unknown_quantities()[i], A, b);

              VectorType x = solve(A, b);

              if (this->quantities().unknown_quantities()[i].get_name() == viennashe::quantity::potential())
              {
                potential_norm_increment = viennashe::math::norm_2(x);
              }
              current_residual_norm += viennashe::math::norm_2(b);

              total_update_norm += viennashe::get_relative_update_norm(this->quantities().unknown_quantities()[i], x);

              viennashe::update_quantity(this->device(), this->quantities().unknown_quantities()[i], this->config(), x, update_no_damping);
            }

            // SHE quantities:
            for (std::size_t i = 0; i < quantities().unknown_she_quantities().size(); ++i)
            {
              typename map_info_type::value_type::second_type quan_unknowns = map_info[this->quantities().unknown_she_quantities()[i].get_name()];
              std::size_t number_of_unknowns = quan_unknowns.first + quan_unknowns.second;
              if (number_of_unknowns == 0)
                continue;

              // System for this quantity only:
              MatrixType A(number_of_unknowns, number_of_unknowns);
              VectorType b(number_of_unknowns);

              viennashe::she::assemble(device(), transferred_quantities, this->quantities(), this->config(),
                                        this->quantities().unknown_she_quantities()[i], A, b,
                                        (quantities_history_.size() > 1), nonlinear_iter > 1);

              VectorType x = viennashe::she::solve(this->device(), this->quantities().unknown_she_quantities()[i], this->config(), A, b, quan_unknowns.first);

              current_residual_norm += viennashe::math::norm_2(b);


              // TODO: The following is a rather ugly dispatch:
              if (quantities().unknown_she_quantities()[i].get_name() == viennashe::quantity::electron_distribution_function())
              {
                total_update_norm += viennashe::get_relative_update_norm(device(), quantities().unknown_she_quantities()[i],
                                    quantities().get_unknown_quantity(viennashe::quantity::electron_density()), config());

                viennashe::update_quantity(device(), quantities().unknown_she_quantities()[i],
                                           quantities().get_unknown_quantity(viennashe::quantity::electron_density()),
                                           config(), x, update_no_damping);
              }
              else
              {
                total_update_norm += viennashe::get_relative_update_norm(device(), quantities().unknown_she_quantities()[i],
                                      quantities().get_unknown_quantity(viennashe::quantity::hole_density()), config());

                viennashe::update_quantity(device(), quantities().unknown_she_quantities()[i],
                                           quantities().get_unknown_quantity(viennashe::quantity::hole_density()),
                                           config(), x, update_no_damping);
              }
            } // for each SHE quantity
          }

          // Print norms:
          log::info<log_simulator>() << std::scientific << std::setprecision(3) << std::setw(9) << current_residual_norm    << " | "
                                     << std::scientific << std::setprecision(3) << std::setw(9) << potential_norm_increment << " | "
                                     << std::scientific << std::setprecision(3) << std::setw(9) << total_update_norm        << " | "
                                     << std::fixed      << std::setprecision(3) << std::setw(8) << stopwatch.get() << std::endl;

          // push every nonlinear iteration state:
          //quantities_history_.push_back(quantities());

          if (break_after_update)
            break;

          if (!current_residual_norm) //equilibrium case: no iterations required
            break;

          // Check for convergence:
          if (nonlinear_iter == 1)
          {
            initial_residual_norm = current_residual_norm;
          }
          else if (current_residual_norm / initial_residual_norm < 1e-10)
          {
            update_no_damping = true;
            break_after_update = true;
            //break;
          }
          else if (potential_norm_increment < this->config().nonlinear_solver().tolerance())
          {
            update_no_damping = true;
            break_after_update = true;
            //break;
          }

        } //for nonlinear iterations


      }


      /** @brief Returns the controller object. Const version. */
      SHETimeStepQuantitiesT const & quantities() const { return quantities_history_.back(); }

      /** @brief Returns the controller object. Non-const version. */
      SHETimeStepQuantitiesT & quantities() { return quantities_history_.back(); }

      /** @brief Returns a wrapper for the distribution function for electrons and holes (whatever is computed by the simulator),
       *         which can be evaluated in the vertices and edges of the mesh */
      she_df_type she_df(carrier_type_id ctype) const
      {
        return she_df_type(config(), this->quantities().carrier_distribution_function(ctype));
      }

      /** @brief Returns a wrapper for the evaluation of the energy distribution function (i.e. the isotropic part of f) */
      edf_type edf(carrier_type_id ctype) const
      {
        return edf_type(config(), this->quantities().carrier_distribution_function(ctype));
      }

      /** @brief Returns a wrapper for the evaluation of the generalized energy distribution function (f * Z, with density of states Z) */
      generalized_edf_type generalized_edf(carrier_type_id ctype) const
      {
        return generalized_edf_type(config(), this->quantities().carrier_distribution_function(ctype));
      }

      /** @brief Returns a wrapper for the potential, which can be evaluated in every vertex of the mesh
       *
       * A simulation must have been carried out already, otherwise a quantity_not_yet_available_exception is thrown.
       */
      potential_type potential() const
      {
        return quantities().potential();
      }

      /** @brief Returns a wrapper for the electron density, which can be evaluated in every vertex of the mesh */
      electron_density_type electron_density() const
      {
        return quantities().electron_density();
      }

      /** @brief Returns a wrapper for the hole density, which can be evaluated in every vertex of the mesh   */
      hole_density_type hole_density() const
      {
        return quantities().hole_density();
      }

      /** @brief Returns a wrapper for the trap occupancy, which can be evaluated in every vertex of the mesh */
      //trap_occupancy_type trap_occupancy() const
      //{
      //  return quantities().trap_occupancy();
      //}

      /** @brief Returns the state of the simulator at an earlier timestep.
       *
       *  @param index  Distance from first time step. 0: first timestep, 1: next timestep, etc.
       */
      SHETimeStepQuantitiesT const & quantity_history(std::size_t index) const { return quantities_history_.at(index); }

      std::size_t quantity_history_size() const { return quantities_history_.size(); }

      /** @brief Cycles the quantities. Moves (copy) the current quantities to the history and empties the current quantities. Does not increment time.   */
      void advance_in_time()
      {
        quantities_history_.push_back(quantities());
        detail::set_boundary_for_material(device(), quantities().get_unknown_quantity(viennashe::quantity::potential()), materials::checker(MATERIAL_CONDUCTOR_ID), boundary_potential_accessor<DeviceType>(device()), BOUNDARY_DIRICHLET);
      }

      // the following is for compatibility reasons:
      ResultQuantityType dg_pot_n() const { return quantities().dg_pot_n(); }
      ResultQuantityType dg_pot_p() const { return quantities().dg_pot_p(); }

      /** @brief Returns the config object used by the simulator controller */
      viennashe::config const & config() const { return config_; }
      viennashe::config       & config()       { return config_; }

      DeviceType const & device() const { return *p_device_; }
      DeviceType & device() { return *p_device_; }

    private:

      /** @brief Transfers the provided source quantities to the destination for the provided elements. */
      template <typename QuantitySrcT, typename QuantityDestT>
      void transfer_provided_quantities(QuantitySrcT const & src_quantity,
                                        QuantityDestT      & dest_quantity)
      {
        // Usage example: Transfer initial guess to SHE index space

        viennagrid_dimension cell_dim;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(p_device_->mesh(), &cell_dim));

        viennagrid_element_id *cells_begin, *cells_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(p_device_->mesh(), cell_dim, &cells_begin, &cells_end));
        for (viennagrid_element_id *cit  = cells_begin;
                                    cit != cells_end;
                                  ++cit)
        {
          dest_quantity.set_value(*cit, src_quantity.at(*cit));
        }

      } // transfer_provided_quantities


      /** @brief Solves the provided linear system by first scaling the unknowns appropriately (preconditioning purposes)
       * @tparam MatrixType The matrix type (cf. linalg_util.hpp).
       * @param matrix The system matrix (sparse matrix)
       * @param rhs The right hand side (usually a std::vector)
       * @param rescaling_vector The scaling vector. Can be empty if scaling is disabled
       */
      template <typename MatrixType>
      VectorType solve(MatrixType & matrix, VectorType & rhs, VectorType & rescaling_vector)
      {
        if (this->config_.linear_solver().scale())
        {
          // Step 1: Rescale unknowns:
          typedef typename MatrixType::row_type      RowType;
          typedef typename MatrixType::iterator2     AlongRowIterator;

          for (std::size_t i=0; i<matrix.size1(); ++i)
          {
            RowType & row_i = matrix.row(i);
            for (AlongRowIterator iter  = row_i.begin();
                                  iter != row_i.end();
                                ++iter)
            {
              iter->second *= rescaling_vector[iter->first];
            }
          }
        }

        // Step 2: Normalize matrix rows:
        viennashe::math::row_normalize_system(matrix, rhs);

        // Step 3: Solve
        VectorType update = viennashe::solvers::solve(matrix, rhs, config().linear_solver());

        // Step 4: Check convergence:
        double lin_sol_res = viennashe::math::norm_2(viennashe::math::subtract(viennashe::math::prod(matrix, update), rhs));
        lin_sol_res       /= viennashe::math::norm_2(rhs);

        if (lin_sol_res > 1e-3)
          log::warning() << "Warning: Linear solver shows only mild convergence! Residual: " << lin_sol_res << std::endl;
        if (lin_sol_res > 1)
        {
          log::error() << "ERROR: Linear solver failed to converge properly! Residual: " << lin_sol_res << std::endl;
          //log::debug() << rhs << std::endl;
          throw viennashe::solver_failed_exception("Linear solver failed to converge properly!");
        }

        if (this->config_.linear_solver().scale())
        {
          // Step 4: Revert unknown scaling:
          for (std::size_t i=0; i<update.size(); ++i)
            update[i] *= rescaling_vector[i];
        }

        return update;
      } // solve

      /**
       *
       * @param matrix The system matrix A
       * @param rhs The right hand side b
       * @return The solution x of A x = b
       */
      template <typename MatrixType>
      VectorType solve(MatrixType & matrix, VectorType & rhs, bool is_newton = false)
      {
        VectorType rescaling_vector(rhs.size());
        if (this->config_.linear_solver().scale())
        {
          std::fill(rescaling_vector.begin(), rescaling_vector.end(), 1.0);

          if (is_newton)
          {

            viennagrid_mesh mesh = device().mesh();

            viennagrid_dimension cell_dim;
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh, &cell_dim));

            viennagrid_element_id *cells_begin, *cells_end;
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end));

            for (std::size_t i=0; i<quantities().unknown_quantities().size(); ++i)
            {
              UnknownQuantityType & current_quan = quantities().unknown_quantities()[i];

              double scaling_factor = 1.0;
              if (   current_quan.get_name() == viennashe::quantity::electron_density()
                  || current_quan.get_name() == viennashe::quantity::hole_density())
                scaling_factor = 1e16;


              for (viennagrid_element_id *cit  = cells_begin;
                                          cit != cells_end;
                                        ++cit)
              {
                long idx = current_quan.get_unknown_index(*cit);
                if (idx >= 0)
                  rescaling_vector[static_cast<std::size_t>(idx)] = scaling_factor;
              }
            }
          }
        }

        return solve(matrix, rhs, rescaling_vector);
      }

      // Transferred from solver.hpp
      /*
      double compute_newton_damping_factor(DeviceType const & device, VectorType const & she_update, VectorType const & n_old, VectorType const & p_old)
      {
        //VectorType
        VectorType she_result_new(she_result_);
        for (std::size_t i=0; i<she_result_.size(); ++i)
          she_result_new[i] += she_update[i];

        VectorType n_new(n_old);
        if (conf_.with_electrons())
          compute_carrier_density_vector(device, quantities(), quantities().electron_distribution_function(), conf_.dispersion_relation(viennashe::ELECTRON_TYPE_ID), n_new);

        VectorType p_new(p_old);
        if (conf_.with_holes())
          compute_carrier_density_vector(device, quantities(), quantities().hole_distribution_function(), conf_.dispersion_relation(viennashe::HOLE_TYPE_ID), p_new);


        // compute relative norm of update:
        double norm_inf_update = 0;
        double norm_inf_current = 0;

        for (std::size_t i=0; i<n_new.size(); ++i)
        {
            norm_inf_update  = std::max(std::abs(norm_inf_update),  std::abs(n_new[i]));
            norm_inf_current = std::max(std::abs(norm_inf_current), std::abs(n_old[i]));
        }

        double norm_rel_update = norm_inf_update / norm_inf_current;

        double delta = 0.1;
        return (1.0 + delta * std::log(norm_rel_update)) / (1.0 + delta * (norm_rel_update - 1.0));


        std::size_t damping_iter_max = 10;
        double alpha = 1.0;

        for (std::size_t damping_iter = 0; damping_iter < damping_iter_max; ++damping_iter)
        {
          for (std::size_t i=0; i<she_result_.size(); ++i)
            she_result_new[i] = she_result_[i] + alpha * she_update[i];

          if (damping_iter == 0) //check whether a full Newton step is okay
          {
            alpha = 1.0;
          }
          else if (damping_iter == 1) //use logarithmic damping, cf. PhD thesis by Simlinger or Fischer
          {
            alpha = (1.0 + delta * std::log(norm_rel_update)) / (1.0 + delta * (norm_rel_update - 1.0));
          }
          else
          {
            alpha = 0.3 * alpha;
          }

          double update_norm = 0.0;
          double current_norm = 0.0;

          bool do_reject = false;

          if (conf_.with_electrons())
            compute_carrier_density_vector(device, quantities(), quantities().electron_distribution_function(), conf_.dispersion_relation(viennashe::ELECTRON_TYPE_ID), n_new);

          if (conf_.with_holes())
            compute_carrier_density_vector(device, quantities(), quantities().hole_distribution_function(), conf_.dispersion_relation(viennashe::HOLE_TYPE_ID), p_new);

          VectorType n_update = n_old;
          for (std::size_t i=0; i<n_update.size(); ++i)
            n_update[i] -= n_new[i];


          for (std::size_t i=0; i<n_new.size(); ++i)
          {
            update_norm += alpha * alpha * n_update[i] * n_update[i];
            current_norm += n_old[i] * n_old[i];

            if (   std::abs(n_old[i]) != 0.0
                && std::abs(n_old[i]) < std::abs(alpha * n_update[i])
                )
              do_reject = true;

          }

          if (update_norm > current_norm)
            do_reject = true;

          // Don't accept updated variables if update is too large
          if ( (damping_iter < damping_iter_max - 1) //last iterate is always accepted
              && do_reject)
          {
            continue;
          }

          break;
        } //for damping_iter

        return alpha;
      } */



      /** @brief Generic helper function for updating a quantity (usually potential or densities). Currently used within the Newton scheme. */
      /*template <typename IndexAccessor, typename DamperType>
      VectorType compute_next_iterate(DeviceType const & device,
                                      VectorType const & update, VectorType const & current_guess, IndexAccessor const & index_accessor,
                                      DamperType const & damper, bool strict_positivity_required)
      {
        //typedef typename DeviceType::mesh_type                                          MeshType;
        typedef typename viennagrid::result_of::const_ncell_range<MeshType, 0>::type      VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type           VertexIterator;

        MeshType const & mesh = device.mesh();

        VectorType next_iterate(current_guess.size());

        // compute relative norm of update:
        double norm_inf_update = 0;
        double norm_inf_current = 0;

        VertexContainer vertices = viennagrid::ncells<0>(mesh);
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {
          long index = index_accessor(*vit);
          if (index >= 0)
          {
            norm_inf_update  = std::max(std::abs(norm_inf_update),  std::abs(update[index]));
            norm_inf_current = std::max(std::abs(norm_inf_current), std::abs(current_guess[index]));
          }
        }

        double norm_rel_update = norm_inf_update / norm_inf_current;

        double delta = 0.1;

        std::size_t damping_iter_max = 10;
        double alpha = 1.0;

        for (std::size_t damping_iter = 0; damping_iter < damping_iter_max; ++damping_iter)
        {
          next_iterate.clear();
          next_iterate.resize(pot_n_p_.size());

          if (damping_iter == 0) //check whether a full Newton step is okay
          {
            alpha = 1.0;
          }
          else if (damping_iter == 1) //use logarithmic damping, cf. PhD thesis by Simlinger or Fischer
          {
            if ( viennashe::util::is_Inf(norm_rel_update) ||
                 viennashe::util::is_NaN(norm_rel_update)) { alpha = damper(alpha); } // limes norm_rel_update -> inf = 1.0
            else { alpha = (1.0 + delta * std::log(norm_rel_update)) / (1.0 + delta * (norm_rel_update - 1.0)); }
          }
          else
          {
            alpha = damper(alpha);
          }

          double update_norm = 0.0;
          double current_norm = 0.0;

          bool do_reject = false;

          // Compute updated potential, electron density and hole density
          VertexContainer vertices = viennagrid::ncells<0>(mesh);
          for (VertexIterator vit = vertices.begin();
              vit != vertices.end();
              ++vit)
          {
            long index = index_accessor(*vit);
            if (index > -1)
            {
              next_iterate[index] = current_guess[index] + alpha * update[index];
              update_norm += alpha * alpha * update[index] * update[index];
              current_norm += current_guess[index] * current_guess[index];

              if (   std::abs(current_guess[index]) != 0.0
                  && std::abs(current_guess[index]) < 0.5 * std::abs(alpha * update[index])
                 )
                do_reject = true;

              if (strict_positivity_required && next_iterate[index] <= 0)
                do_reject = true;
            }
          }

          if (update_norm > current_norm)
            do_reject = true;

          // Don't accept updated variables if update is too large
          if ( (damping_iter < damping_iter_max - 1) //last iterate is always accepted
              && do_reject)
          {
            continue;
          }

          break;
        } //for update

        return next_iterate;
      } */

      DeviceType * p_device_;
      viennashe::config config_;

      std::vector<SHETimeStepQuantitiesT> quantities_history_;

  }; //simulator

} //namespace viennashe

#endif
