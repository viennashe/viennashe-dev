#ifndef VIENNASHE_SHE_TIMESTEP_QUANTITIES_HPP
#define VIENNASHE_SHE_TIMESTEP_QUANTITIES_HPP

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

// viennashe
#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/config.hpp"
#include "viennashe/she/df_wrappers.hpp"
#include "viennashe/she/she_quantity.hpp"
#include "viennashe/accessors.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/exception.hpp"

/** @file viennashe/she/timestep_quantities.hpp
    @brief A container of all quantities defined for a certain timestep t.
*/

namespace viennashe
{
  namespace she
  {


    /** @brief The main SHE simulator controller class. Acts as an accessor for all SHE quantities needed for the simulation. */
    template <typename DeviceType>
    class timestep_quantities
    {
        typedef timestep_quantities<DeviceType>     self_type;
        typedef std::vector<double>                 VectorType;

        typedef std::vector<std::vector<long> >    index_vector_type;
        typedef std::vector<std::vector<double> >  trap_occupancy_type;

        typedef std::vector<long>                 expansion_order_vector_type;
        typedef std::vector<double>               boundary_values_vector_type;
        typedef std::vector<double>               error_values_vector_type;

      public:
        typedef DeviceType               device_type;

        typedef unknown_she_quantity<viennagrid_element_id, viennagrid_element_id>     UnknownSHEQuantityType;
        typedef UnknownSHEQuantityType                                                 unknown_she_quantity_type;
        typedef std::deque<UnknownSHEQuantityType>                                     UnknownSHEQuantityListType;
        typedef const_she_quantity<viennagrid_element_id, viennagrid_element_id>       ResultSHEQuantityType;

        typedef viennashe::she::she_df_wrapper<DeviceType, UnknownSHEQuantityType>                       she_df_type;
        typedef viennashe::she::edf_wrapper<DeviceType, UnknownSHEQuantityType>                             edf_type;
        typedef viennashe::she::generalized_edf_wrapper<DeviceType, UnknownSHEQuantityType>     generalized_edf_type;

        typedef unknown_quantity<viennagrid_element_id>       UnknownQuantityType;
        typedef UnknownQuantityType                           unknown_quantity_type;
        typedef const_quantity<viennagrid_element_id>         ResultQuantityType;

        typedef ResultQuantityType             potential_type;
        typedef ResultQuantityType             electron_density_type;
        typedef ResultQuantityType             hole_density_type;


        // use default implementations (should suffice)
        //timestep_quantities() {}
        //simulator_controller(simulator_controller const &) {}
        //void operator=(simulator_controller const &) {}


        /** @brief Returns the number of inner unknown indices (aka. degree of freedom - dof) for traps associated with a cell.
         *
         * @param c          The respective cell
         */
        std::size_t num_trap_unknown_indices(viennagrid_element_id c) const
        {
          return cell_trap_unknown_indices_.at(std::size_t(viennagrid_index_from_element_id(c))).size();
        }

        /** @brief Returns an unknown index for traps associated with a cell
         *
         * @param c            The respective cell
         * @param inner_index  The inner index
         */
        long trap_unknown_index(viennagrid_element_id c, long inner_index) const
        {
          return cell_trap_unknown_indices_.at(viennagrid_index_from_element_id(c)).at(inner_index);
        }

        void setup_trap_unkown_indices(DeviceType const & device)
        {
          typedef typename DeviceType::trap_level_container_type     TrapContainerType;
          typedef typename TrapContainerType::const_iterator         TrapIterator;

          long i = 0;
          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

          viennagrid_element_id *cells_begin, *cells_end;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end));

          cell_trap_unknown_indices_.resize(cells_end - cells_begin);
          cell_trap_occupancies_.resize(cells_end - cells_begin);

          for (viennagrid_element_id *cit  = cells_begin;
                                      cit != cells_end;
                                    ++cit)
          {
            std::size_t cell_index(viennagrid_index_from_element_id(*cit));
            TrapContainerType const & traps = device.get_trap_levels(*cit);

            if (viennashe::materials::is_semiconductor(device.get_material(*cit)) && traps.size() > 0 )
            {
              std::vector<long> & idx = cell_trap_unknown_indices_.at(cell_index);
              idx.resize(traps.size());
              cell_trap_occupancies_.at(cell_index).resize(traps.size(), 0);

              std::size_t j = 0;
              for (TrapIterator trap_it  = traps.begin(); trap_it != traps.end(); ++trap_it, ++i, ++j)
                idx[j] = i;
            }
            else
            {
              // NO UNKNOWNS HERE !
              cell_trap_unknown_indices_.at(cell_index).clear();
              cell_trap_unknown_indices_.at(cell_index).resize(0);
            }
          }
        }

        double trap_occupancy(viennagrid_element_id c, std::size_t inner_index) const
        {
          return cell_trap_occupancies_.at(std::size_t(viennagrid_index_from_element_id(c))).at(inner_index);
        }

        void trap_occupancy(viennagrid_element_id c, std::size_t inner_index, double new_occupancy)
        {
          if(new_occupancy < 0.0 || new_occupancy > 1.0)
          {
            log::error() << "ERROR: Invalid trap occupancy: " << new_occupancy << std::endl;
            throw viennashe::invalid_value_exception("trap_level.occupancy: occupancies have to be between 0 and 1!", new_occupancy);
          }
          cell_trap_occupancies_.at(std::size_t(viennagrid_index_from_element_id(c))).at(inner_index) = new_occupancy;
        }

        ////////////// SHE quantities ////////////////////////

        /** @brief Returns a const reference to a SHE quantity identified by its name.
         *
         * @param quantity_name A std::string uniquely identifying the quantity
         * @throw quantity_not_found_exception This method may throw a quantity_not_found_exception in case the requested quantity could not be found
         */
        UnknownSHEQuantityType const & she_quantity(std::string quantity_name) const
        {
          for (std::size_t i=0; i<unknown_she_quantities_.size(); ++i)
          {
            if (quantity_name == unknown_she_quantities_[i].get_name())
              return unknown_she_quantities_[i];
          }

          // quantity not found -> throw exception
          std::stringstream ss;
          ss << "No SHE quantity named '" << quantity_name << "'";
          throw quantity_not_found_exception(ss.str());
        }

        /** @brief Returns a reference to a SHE quantity identified by its name.
         *
         * @param quantity_name A std::string uniquely identifying the quantity
         * @throw quantity_not_found_exception This method may throw a quantity_not_found_exception in case the requested quantity could not be found
         */
        UnknownSHEQuantityType & she_quantity(std::string quantity_name)
        {
          for (std::size_t i=0; i<unknown_she_quantities_.size(); ++i)
          {
            if (quantity_name == unknown_she_quantities_[i].get_name())
              return unknown_she_quantities_[i];
          }

          // quantity not found -> throw exception
          std::stringstream ss;
          ss << "No SHE quantity named '" << quantity_name << "'";
          throw quantity_not_found_exception(ss.str());
        }

        UnknownSHEQuantityType const &  electron_distribution_function() const { return she_quantity(viennashe::quantity::electron_distribution_function()); }
        UnknownSHEQuantityType       &  electron_distribution_function()       { return she_quantity(viennashe::quantity::electron_distribution_function()); }

        UnknownSHEQuantityType const &      hole_distribution_function() const { return she_quantity(viennashe::quantity::hole_distribution_function()); }
        UnknownSHEQuantityType       &      hole_distribution_function()       { return she_quantity(viennashe::quantity::hole_distribution_function()); }

        UnknownSHEQuantityType const &      carrier_distribution_function(viennashe::carrier_type_id ctype) const
        {
          if (ctype == viennashe::ELECTRON_TYPE_ID)  return this->electron_distribution_function();
          else if (ctype == viennashe::HOLE_TYPE_ID) return this->hole_distribution_function();
          else throw viennashe::carrier_type_not_supported_exception("In carrier_distribution_function");
        }
        UnknownSHEQuantityType       &      carrier_distribution_function(viennashe::carrier_type_id ctype)
        {
          return static_cast<const self_type*>(this)->carrier_distribution_function(ctype); // call the const version; avoid code duplication
        }

        UnknownSHEQuantityListType       & unknown_she_quantities()       { return unknown_she_quantities_; }
        UnknownSHEQuantityListType const & unknown_she_quantities() const { return unknown_she_quantities_; }

        ////////////// Macroscopic quantities ////////////////////////

        /** @brief Returns the quantity identified by its name.
         *
         * @param quantity_name A std::string uniquely identifying the quantity
         * @throw quantity_not_found_exception This method may throw a quantity_not_found_exception in case the requested quantity could not be found
         */
        ResultQuantityType quantity(std::string quantity_name) const
        {
          for (std::size_t i=0; i<unknown_quantities_.size(); ++i)
          {
            if (quantity_name == unknown_quantities_[i].get_name())
              return ResultQuantityType(quantity_name, unknown_quantities_[i].values());
          }

          // quantity not found -> throw exception
          std::stringstream ss;
          ss << "No quantity named '" << quantity_name << "'";
          throw quantity_not_found_exception(ss.str());
        }

        ResultQuantityType        potential()    const { return quantity(viennashe::quantity::potential()); }
        ResultQuantityType electron_density()    const { return quantity(viennashe::quantity::electron_density()); }
        ResultQuantityType     hole_density()    const { return quantity(viennashe::quantity::hole_density()); }
        ResultQuantityType         dg_pot_n()    const { return quantity(viennashe::quantity::density_gradient_electron_correction()); }
        ResultQuantityType         dg_pot_p()    const { return quantity(viennashe::quantity::density_gradient_hole_correction()); }
        ResultQuantityType lattice_temperature() const { return quantity(viennashe::quantity::lattice_temperature()); }

        std::deque<UnknownQuantityType>       & unknown_quantities()       { return unknown_quantities_; }
        std::deque<UnknownQuantityType> const & unknown_quantities() const { return unknown_quantities_; }

        /** @brief Returns a reference to the <b>unkown</b> quantity identified by its name.
         *
         * @param quantity_name A std::string uniquely identifying the quantity
         * @throw quantity_not_found_exception This method may throw a quantity_not_found_exception in case the requested quantity could not be found
         */
        UnknownQuantityType & get_unknown_quantity(std::string quantity_name)
        {
          for (std::size_t i=0; i<unknown_quantities_.size(); ++i)
          {
            if (quantity_name == unknown_quantities_[i].get_name())
              return unknown_quantities_[i];
          }

          // quantity not found -> throw exception
          std::stringstream ss;
          ss << "No unknown quantity named '" << quantity_name << "'";
          throw quantity_not_found_exception(ss.str());
        }

        /** @brief Returns a const reference to the <b>unkown</b> quantity identified by its name.
         *
         * @param quantity_name A std::string uniquely identifying the quantity
         * @throw quantity_not_found_exception This method may throw a quantity_not_found_exception in case the requested quantity could not be found
         */
        UnknownQuantityType const & get_unknown_quantity(std::string quantity_name) const
        {
          for (std::size_t i=0; i<unknown_quantities_.size(); ++i)
          {
            if (quantity_name == unknown_quantities_[i].get_name())
              return unknown_quantities_[i];
          }

          // quantity not found -> throw exception
          std::stringstream ss;
          ss << "No unknown quantity named '" << quantity_name << "'";
          throw quantity_not_found_exception(ss.str());
        }

      private:

        // Quantities in SHE space (at least f^n and/or f^p):
        UnknownSHEQuantityListType      unknown_she_quantities_;

        // Poisson and DD related:
        std::deque<UnknownQuantityType> unknown_quantities_;

        // traps ...
        index_vector_type               cell_trap_unknown_indices_;
        trap_occupancy_type             cell_trap_occupancies_;


    };

  }
}

#endif
