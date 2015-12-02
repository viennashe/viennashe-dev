#ifndef VIENNASHE_SRH_KINETICS_HPP
#define	VIENNASHE_SRH_KINETICS_HPP
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
#include "viennashe/log/log.hpp"

#include "viennashe/trap_level.hpp"
#include "viennashe/she/assemble_common.hpp"

namespace viennashe
{
  namespace models
  {

    namespace srh
    {

      namespace detail
      {

        template < typename DeviceType, typename TimeStepQuantitiesT >
        double get_carrier_concentration(DeviceType const & device,
                                         viennagrid_element_id cell_or_facet,
                                         TimeStepQuantitiesT const & quantities,
                                         viennashe::carrier_type_id ctype)
        {
          typedef typename TimeStepQuantitiesT::ResultQuantityType  CarrierDensityType;

          viennagrid_dimension cell_dim;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim));

          CarrierDensityType const & carrier_density = (ctype == viennashe::ELECTRON_TYPE_ID) ? quantities.electron_density() : quantities.hole_density();

          if (cell_dim == viennagrid_topological_dimension_from_element_id(cell_or_facet)) // cell
          {
              return carrier_density.at(cell_or_facet);
          }
          else if (cell_dim == viennagrid_topological_dimension_from_element_id(cell_or_facet) + 1) // facet
          {
            viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
            VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(device.mesh(), cell_or_facet, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

            if (cells_on_facet_begin + 1 == cells_on_facet_end)
              return carrier_density.at(cells_on_facet_begin[0]);

            if (!viennashe::materials::is_semiconductor(device.get_material(cells_on_facet_begin[0])))
              return carrier_density.at(cells_on_facet_begin[1]);
            else if (!viennashe::materials::is_semiconductor(device.get_material(cells_on_facet_begin[1])))
              return carrier_density.at(cells_on_facet_begin[0]);
            else
            {
              viennagrid_numeric density_0 = carrier_density.at(cells_on_facet_begin[0]);
              viennagrid_numeric density_1 = carrier_density.at(cells_on_facet_begin[1]);
              return std::sqrt(density_0 * density_1);
            }
          }
          else
            throw std::runtime_error("get_carrier_concentration(): Invalid element type!");
        }


        /**
         * @brief Implementation of the recombination term without any occupancies considered
         * @param trap The SRH trap level
         * @param device The device
         * @param el The element for which to calculate the recombination
         * @param conf The simulator configuration
         * @param quantities The timestep quantities
         * @param ctype The carrier type
         * @param index_H If SHE is being used: the H-space index
         * @return A recombination term, either just spatial (DD) or per energy (SHE)
         */
        template < typename DeviceType, typename ElementType, typename TimeStepQuantitiesT >
        double gamma_recombination_impl(viennashe::trap_level const & trap,
                                        DeviceType const & device,
                                        ElementType const & el,
                                        viennashe::config const & conf,
                                        TimeStepQuantitiesT const & quantities,
                                        viennashe::carrier_type_id ctype,
                                        std::size_t index_H = 0)
        {
          const double collision_cs = trap.collision_cross_section();

          if (ctype == viennashe::ELECTRON_TYPE_ID && conf.with_electrons() && conf.get_electron_equation() == viennashe::EQUATION_SHE)
          {
            if (conf.get_electron_equation() == viennashe::EQUATION_SHE)
            {
              for (std::size_t i = 0; i < quantities.unknown_she_quantities().size(); ++i)
              {
               // TODO: The following is a rather ugly dispatch:
               if (quantities.unknown_she_quantities()[i].get_name() != viennashe::quantity::electron_distribution_function())
                 continue;

                const double kinetic_energy = quantities.unknown_she_quantities()[i].get_kinetic_energy(el, index_H);
                const double velocity       = conf.dispersion_relation(ctype).velocity(kinetic_energy);
                double gamma_recombination  = collision_cs * velocity;
                return gamma_recombination;
              }
            }
            else if (conf.get_electron_equation() == viennashe::EQUATION_CONTINUITY)
            {
              const double  T  = device.get_lattice_temperature(el);
              const double vth = viennashe::physics::get_thermal_velocity(T, viennashe::ELECTRON_TYPE_ID);
              const double sigma_n = vth * collision_cs;

              return sigma_n * get_carrier_concentration(device, el, quantities, viennashe::ELECTRON_TYPE_ID);
            }
          }
          else if (ctype == viennashe::HOLE_TYPE_ID && conf.with_holes())
          {
            if ( conf.get_hole_equation() == viennashe::EQUATION_SHE )
            {
              for (std::size_t i = 0; i < quantities.unknown_she_quantities().size(); ++i)
              {
               // TODO: The following is a rather ugly dispatch:
               if (quantities.unknown_she_quantities()[i].get_name() != viennashe::quantity::hole_distribution_function())
                 continue;

                const double kinetic_energy = quantities.unknown_she_quantities()[i].get_kinetic_energy(el, index_H);
                const double velocity       = conf.dispersion_relation(ctype).velocity(kinetic_energy);
                double gamma_recombination  = collision_cs * velocity;

                return gamma_recombination;
              }
            }
            else if (conf.get_hole_equation() == viennashe::EQUATION_CONTINUITY)
            {
              const double  T  = device.get_lattice_temperature(el);
              const double vth = viennashe::physics::get_thermal_velocity(T, viennashe::HOLE_TYPE_ID);
              const double sigma_p = vth * collision_cs;

              return sigma_p * get_carrier_concentration(device, el, quantities, viennashe::HOLE_TYPE_ID);
            }
          }

          throw viennashe::unavailable_feature_exception("SRH-trap evaluation needs the electron and hole EDF or electron and hole concentrations!");

          return 0;
        }

        /**
         * @brief Implementation of the generation term without any occupancies considered
         * @param trap The SRH trap level
         * @param device The device
         * @param el The element for which to calculate the generation
         * @param conf The simulator configuration
         * @param quantities The timestep quantities
         * @param ctype The carrier type
         * @param index_H If SHE is being used: the H-space index
         * @return A generation term, either just spatial (DD) or per energy (SHE)
         */
        template < typename DeviceType, typename ElementType, typename TimeStepQuantitiesT >
        double gamma_generation_impl(viennashe::trap_level const & trap,
                                     DeviceType const & device,
                                     ElementType const & el,
                                     viennashe::config const & conf,
                                     TimeStepQuantitiesT const & quantities,
                                     viennashe::carrier_type_id ctype,
                                     std::size_t index_H = 0)
        {
          const double collision_cs     = trap.collision_cross_section();
          const double cell_temperature = device.get_lattice_temperature(el);
          const double kBT              = viennashe::physics::constants::kB * cell_temperature;

          if (ctype == viennashe::ELECTRON_TYPE_ID && conf.with_electrons())
          {
            if (conf.get_electron_equation() == viennashe::EQUATION_SHE)
            {
              for (std::size_t i = 0; i < quantities.unknown_she_quantities().size(); ++i)
              {
                // TODO: The following is a rather ugly dispatch:
                if (quantities.unknown_she_quantities()[i].get_name() != viennashe::quantity::electron_distribution_function())
                  continue;

                double total_trap_energy    = trap.energy() + quantities.unknown_she_quantities()[i].get_bandedge_shift(el) - viennashe::physics::get_band_edge(viennashe::ELECTRON_TYPE_ID); //total trap energy
                const double kinetic_energy = quantities.unknown_she_quantities()[i].get_kinetic_energy(el, index_H);
                const double velocity       = conf.dispersion_relation(ctype).velocity(kinetic_energy);
                const double total_energy   = quantities.unknown_she_quantities()[i].get_value_H(index_H);
                double gamma_generation     = collision_cs * velocity * std::exp( (total_trap_energy - total_energy) / kBT );

                return gamma_generation;
              }
            }
            else if (conf.get_electron_equation() == viennashe::EQUATION_CONTINUITY)
            {
              const double vth = viennashe::physics::get_thermal_velocity(cell_temperature, viennashe::ELECTRON_TYPE_ID);

              const double sigma_n = vth * collision_cs;
              const double n_aux   = viennashe::physics::get_auxilary_concentration(cell_temperature, viennashe::ELECTRON_TYPE_ID);

              return sigma_n * n_aux * std::exp( + trap.energy() / kBT);
            }
          }
          else if (ctype == viennashe::HOLE_TYPE_ID && conf.with_holes())
          {
            if (conf.get_hole_equation() == viennashe::EQUATION_SHE)
            {
              for (std::size_t i = 0; i < quantities.unknown_she_quantities().size(); ++i)
              {
                // TODO: The following is a rather ugly dispatch:
                if (quantities.unknown_she_quantities()[i].get_name() != viennashe::quantity::hole_distribution_function())
                  continue;

                double total_trap_energy    = trap.energy() + quantities.unknown_she_quantities()[i].get_bandedge_shift(el) - viennashe::physics::get_band_edge(viennashe::HOLE_TYPE_ID); //total trap energy
                const double kinetic_energy = quantities.unknown_she_quantities()[i].get_kinetic_energy(el, index_H);
                const double velocity       = conf.dispersion_relation(ctype).velocity(kinetic_energy);
                const double total_energy   = quantities.unknown_she_quantities()[i].get_value_H(index_H);
                double gamma_generation     = collision_cs * velocity * std::exp( (total_energy - total_trap_energy) / kBT ) ;

                return gamma_generation;
              }
            }
            else if (conf.get_hole_equation() == viennashe::EQUATION_CONTINUITY)
            {
              const double vth = viennashe::physics::get_thermal_velocity(cell_temperature, viennashe::HOLE_TYPE_ID);

              const double sigma_p = vth * collision_cs;
              const double p_aux   = viennashe::physics::get_auxilary_concentration(cell_temperature, viennashe::HOLE_TYPE_ID);

              return sigma_p * p_aux * std::exp( - trap.energy() / kBT);
            }
          }

          throw viennashe::unavailable_feature_exception("SRH-trap evaluation needs the electron and hole EDF!");

          return 0;
        }

      } // namespace detail

      /**
       * @brief Returns the trap occupancy based on a bipolar SHE or DD solution
       * @param trap The SRH trap level
       * @param device The device
       * @param el The element on which to operate
       * @param conf The simulator configuration
       * @param quantities The timestep quantities
       */
      template < typename DeviceType, typename ElementType, typename TimeStepQuantitiesT >
      double evaluate(viennashe::trap_level const & trap,
                      DeviceType const & device,
                      ElementType const & el,
                      viennashe::config const & conf,
                      TimeStepQuantitiesT const & quantities)
      {
        const bool with_full_newton = (conf.nonlinear_solver().id() == viennashe::solvers::nonlinear_solver_ids::newton_nonlinear_solver);

        // TODO: We assume SHE in here !

        if (!conf.with_electrons() || conf.get_electron_equation() != viennashe::EQUATION_SHE)
          throw viennashe::quantity_not_found_exception("SRH-trap evaluation needs the electron EDF!");
        if (!conf.with_holes() || conf.get_hole_equation() != viennashe::EQUATION_SHE)
          throw viennashe::quantity_not_found_exception("SRH-trap evaluation needs the hole EDF!");

        // Check for Newton
        if (with_full_newton)
        {
          throw viennashe::unavailable_feature_exception("SRH-trap evaluation not impelemented for Newton!");
        }
        else // GUMMEL SOLVER
        {
          double electron_rec_rate = 0;
          double electron_gen_rate = 0;

          double hole_rec_rate = 0;
          double hole_gen_rate = 0;

          viennagrid_numeric box_volume;
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_volume(device.mesh(), el, &box_volume));

          typedef typename viennashe::config::dispersion_relation_type  DispersionRelation;

          DispersionRelation dispersion_n = conf.dispersion_relation(viennashe::ELECTRON_TYPE_ID);
          DispersionRelation dispersion_p = conf.dispersion_relation(viennashe::HOLE_TYPE_ID);

          for (std::size_t i = 0; i < quantities.unknown_she_quantities().size(); ++i)
          {
            for (std::size_t index_H = 1; index_H < quantities.unknown_she_quantities()[i].get_value_H_size()-1; ++index_H)
            {
              const long f00_index = quantities.unknown_she_quantities()[i].get_unknown_index(el, index_H);
              if (f00_index < 0) //no DOF here
                continue;

              const double energy_spacing     = viennashe::she::box_height(quantities.unknown_she_quantities()[i], el, index_H);
              const bool   assemble_electrons = (quantities.unknown_she_quantities()[i].get_carrier_type_id() == viennashe::ELECTRON_TYPE_ID);
              const double kinetic_energy     = quantities.unknown_she_quantities()[i].get_kinetic_energy(el, index_H);

              const double f00 = quantities.unknown_she_quantities()[i].get_values(el, index_H)[0];

              if (assemble_electrons)
              {
                const double gamma_recombination = detail::gamma_recombination_impl(trap, device, el, conf, quantities, viennashe::ELECTRON_TYPE_ID, index_H);
                const double gamma_generation    = detail::gamma_generation_impl(trap, device, el, conf, quantities, viennashe::ELECTRON_TYPE_ID, index_H);

                electron_rec_rate += box_volume * energy_spacing * gamma_recombination * f00 * dispersion_n.density_of_states(kinetic_energy);
                electron_gen_rate += box_volume * energy_spacing * gamma_generation    * f00 * dispersion_n.density_of_states(kinetic_energy);
              }
              else
              {
                const double gamma_recombination = detail::gamma_recombination_impl(trap, device, el, conf, quantities, viennashe::HOLE_TYPE_ID, index_H);
                const double gamma_generation    = detail::gamma_generation_impl(trap, device, el, conf, quantities, viennashe::HOLE_TYPE_ID, index_H);

                hole_rec_rate     += box_volume * energy_spacing * gamma_recombination * f00 * dispersion_p.density_of_states(kinetic_energy);
                hole_gen_rate     += box_volume * energy_spacing * gamma_generation    * f00 * dispersion_p.density_of_states(kinetic_energy);
              }
            } // for index_H
          }
          // Explicit formula for trap occupancy, cf. SISPAD paper 'Bipolar Spherical Harmonics Expansions of the BTE' by Rupp et al. (2012)
          // Set the trap occupancy
          double ft = 0.0;
          double denominator = electron_rec_rate + electron_gen_rate + hole_rec_rate + hole_gen_rate;
          if (denominator > 0)
            ft = (electron_rec_rate + hole_gen_rate) / denominator;

          // A bit of over- and underflowing occupancy is ok
          if (ft > 1.0 && ft < 1.001)  ft = 1.0;
          if (ft < 0   && ft > -0.001) ft = 0.0;

          //std::cout << el << " => " << ft << std::endl;

          return ft;
        }
      }

      /**
       * @brief Returns the carrier recombination, where the occupancies have been considered!
       * @param trap The SRH trap level
       * @param device The device
       * @param el The element for which to calculate the recombination
       * @param conf The simulator configuration
       * @param quantities The timestep quantities
       * @param ctype The carrier type
       * @param occupancy    Occupancy of the trap (value between 0 and 1)
       * @param index_H If SHE is being used: the H-space index
       * @return A recombination term, either just spatial (DD) or per energy (SHE)
       */
      template < typename DeviceType, typename ElementType, typename TimeStepQuantitiesT >
      double gamma_recombination(viennashe::trap_level const & trap, DeviceType const & device, ElementType const & el, viennashe::config const & conf,
                                 TimeStepQuantitiesT const & quantities,
                                 viennashe::carrier_type_id ctype,
                                 double occupancy,
                                 std::size_t index_H = 0)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return detail::gamma_recombination_impl(trap, device, el, conf, quantities, ctype, index_H) * (1.0 - occupancy);
        else
          return detail::gamma_recombination_impl(trap, device, el, conf, quantities, ctype, index_H) * occupancy;
      }

      /**
       * @brief Returns the carrier generation, where the occupancies have been considered!
       * @param trap The SRH trap level
       * @param device The device
       * @param el The element for which to calculate the generation
       * @param conf The simulator configuration
       * @param quantities The timestep quantities
       * @param ctype The carrier type
       * @param occupancy    Occupancy of the trap (value between 0 and 1)
       * @param index_H If SHE is being used: the H-space index
       * @return A generation term, either just spatial (DD) or per energy (SHE)
       */
      template < typename DeviceType, typename ElementType, typename TimeStepQuantitiesT >
      double gamma_generation(viennashe::trap_level const & trap, DeviceType const & device, ElementType const & el, viennashe::config const & conf,
                              TimeStepQuantitiesT const & quantities,
                              viennashe::carrier_type_id ctype,
                              double occupancy,
                              std::size_t index_H = 0)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return detail::gamma_generation_impl(trap, device, el, conf, quantities, ctype, index_H) * occupancy;
        else
          return detail::gamma_generation_impl(trap, device, el, conf, quantities, ctype, index_H) * (1.0 - occupancy);
      }


    } // namespace models

  } // namespace srh

} // namespace viennashe


#endif	/* VIENNASHE_SRH_KINETICS_HPP */

