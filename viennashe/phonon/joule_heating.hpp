#ifndef VIENNASHE_PHONON_JOULE_HEATING_HPP
#define VIENNASHE_PHONON_JOULE_HEATING_HPP
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

#include "viennashe/forwards.h"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/models/mobility.hpp"
#include "viennashe/postproc/current_density.hpp"
#include "viennashe/postproc/electric_field.hpp"

#include "viennashe/she/postproc/current_density.hpp"

namespace viennashe
{
  /** @brief The namespace for the heat diffusion equation */
  namespace hde
  {
    /**
     * @brief Calculates the Joule Heating (in Watt per Second)
     * @param j The current density accessor
     * @param e The electric field accessor
     */
    template < typename CurrentDensityAccessorType, typename ElectricFieldAccessorType>
    struct joule_heating
    {
       joule_heating(CurrentDensityAccessorType const & j,
                     ElectricFieldAccessorType  const & e )
         : efield(e), jfield(j) { }

       template < typename ElementType >
       double operator()(ElementType const & el) const
       {
         std::vector<double> j = jfield(el);
         std::vector<double> e = efield(el);
         double val = 0.0;
         for (std::size_t i = 0; i < e.size(); ++i)
         {
           val += e[i] * j[i];
         }
         return std::fabs(val);
       }

     private:
       ElectricFieldAccessorType   const & efield;
       CurrentDensityAccessorType  const & jfield;
    };

    /**
     * @brief Auxilary function to joule_heating. Does the same as the joule_heating functor
     * @param jfield An accessor to the current density
     * @param efield An accessor to the electric field
     * @param cell The cell on which the power density is needed
     * @return Power density
     */
    template < typename CurrentDensityAccessorType,
               typename ElectricFieldAccessorType,
               typename CellType >
    double apply_joule_heating(CurrentDensityAccessorType const & jfield,
                             ElectricFieldAccessorType  const & efield,
                             CellType const & cell)
    {
      joule_heating<CurrentDensityAccessorType, ElectricFieldAccessorType> joule_pd(jfield, efield);
      return joule_pd(cell);
    }


   /** @brief Power density accessor. Used to get the power density in the assembly of the heat diffusion equation */
   template <typename DeviceType, typename QuantitiesListType>
   class power_density_accessor
   {
    public:
      typedef typename QuantitiesListType::unknown_quantity_type                quantity_type;
      typedef typename QuantitiesListType::unknown_she_quantity_type            she_quantity_type;

      typedef typename viennashe::models::result_of::mobility_type<DeviceType>::type    mobility_type;

      typedef double    value_type;

      power_density_accessor(DeviceType const & d,
                             QuantitiesListType const & quantities,
                             viennashe::config const & conf)
         : device_(d), conf_(conf),
           mobility_model_n_(viennashe::models::create_constant_mobility_model(d, 0.1430)),
           mobility_model_p_(viennashe::models::create_constant_mobility_model(d, 0.0480)),
           Jfield_dd_n_(d, viennashe::ELECTRON_TYPE_ID, quantities.get_unknown_quantity(viennashe::quantity::potential()), quantities.get_unknown_quantity(viennashe::quantity::electron_density()), mobility_model_n_),
           Jfield_dd_p_(d, viennashe::HOLE_TYPE_ID,     quantities.get_unknown_quantity(viennashe::quantity::potential()), quantities.get_unknown_quantity(viennashe::quantity::hole_density()),         mobility_model_p_),
           Jfield_she_n_(d, conf, quantities.electron_distribution_function()),
           Jfield_she_p_(d, conf, quantities.hole_distribution_function()),
           Efield_(d, quantities.get_unknown_quantity(viennashe::quantity::potential())) {}

      /**
       * @brief Returns the power density depending on the equations for electrons and holes (supports SHE and DD equations)
       * @param cell A cell
       * @return Power density on vertex
       */
      std::vector<value_type> operator()(viennagrid_element_id cell) const
      {
        value_type value = 0;
        // PURE DD
        if ( conf_.get_electron_equation() == viennashe::EQUATION_CONTINUITY &&
             conf_.get_hole_equation()     == viennashe::EQUATION_CONTINUITY )
        {
          value = viennashe::hde::apply_joule_heating(Jfield_dd_n_, Efield_, cell) +
                  viennashe::hde::apply_joule_heating(Jfield_dd_p_, Efield_, cell);
        }
        // PURE SHE (bipolar SHE)
        else if ( conf_.get_electron_equation() == viennashe::EQUATION_SHE &&
                  conf_.get_hole_equation()     == viennashe::EQUATION_SHE )
        {
          value = viennashe::hde::apply_joule_heating(Jfield_she_n_, Efield_, cell) +
                  viennashe::hde::apply_joule_heating(Jfield_she_p_, Efield_, cell);
        }
        // SHE for electrons only
        else if ( conf_.get_electron_equation() == viennashe::EQUATION_SHE &&
                  conf_.get_hole_equation()     == viennashe::EQUATION_CONTINUITY )
        {
          value = viennashe::hde::apply_joule_heating(Jfield_she_n_, Efield_, cell) +
                  viennashe::hde::apply_joule_heating(Jfield_dd_p_,  Efield_, cell);
        }
        else // SHE for holes only
        {
          value = viennashe::hde::apply_joule_heating(Jfield_dd_n_,  Efield_, cell) +
                 viennashe::hde::apply_joule_heating(Jfield_she_p_, Efield_, cell);
        }

        std::vector<value_type> ret(3);
        ret[0] = value;
        return ret;
      }

    private:
      DeviceType         const & device_;
      viennashe::config  const & conf_;

      mobility_type mobility_model_n_;
      mobility_type mobility_model_p_;

      viennashe::current_density_wrapper<DeviceType, quantity_type, quantity_type, mobility_type> Jfield_dd_n_;
      viennashe::current_density_wrapper<DeviceType, quantity_type, quantity_type, mobility_type> Jfield_dd_p_;
      viennashe::she::current_density_wrapper<DeviceType, she_quantity_type>       Jfield_she_n_;
      viennashe::she::current_density_wrapper<DeviceType, she_quantity_type>       Jfield_she_p_;
      viennashe::electric_field_wrapper<DeviceType, quantity_type>                 Efield_;
    };


  } // namespace hde


} // namespace viennashe

#endif // VIENNASHE_PHONON_JOULE_HEATING_HPP

