#ifndef VIENNASHE_SHE_SCATTERING_CONFIG_HPP
#define VIENNASHE_SHE_SCATTERING_CONFIG_HPP
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


#include <iostream>
#include <algorithm>
#include <cmath>

#include "viennashe/forwards.h"
#include "viennashe/log/log.hpp"
#include "viennashe/physics/constants.hpp"

/** @file viennashe/she/scattering/config.hpp
    @brief Scattering operators (and their parameters) are defined here.
*/

namespace viennashe
{

  namespace she
  {

    /** @brief Common base class for all scattering parameter classes. Provides enable/disable interface */
    class scattering_parameter_base
    {
    public:
      scattering_parameter_base(bool b) : enabled_(b) {}

      bool enabled() const { return enabled_; }
      void enabled(bool b) { enabled_ = b; }

    private:
      bool   enabled_;
    };

    /** @brief Parameters for elastic acoustic phonon scattering */
    class acoustic_phonon_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      acoustic_phonon_scattering_parameters() : base_type(true), fit_n_(1.05), fit_p_(0.875) {}

      double get_mass_density(viennashe::carrier_type_id) const { return 2330.0; /* kg/m^3 */ }
      double get_longitudinal_sound_velocity(viennashe::carrier_type_id) const { return 9040.0; /* m/s */  }
      double get_deformation_potential(viennashe::carrier_type_id ctype) const
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return 8.9 * viennashe::physics::constants::q; /* Joule */
        else
          return 5.12 * viennashe::physics::constants::q;
      }

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_   /* vecci: 0.48 */
                                                      : fit_p_;  /* vecci: 0.9 */
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };


    /** @brief Parameters for inelastic optical phonon scattering in single valley approximation */
    class optical_phonon_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      optical_phonon_scattering_parameters() : base_type(true), fit_n_(0.7), fit_p_(0.7) {}

      double get_mass_density(viennashe::carrier_type_id) const { return 2330.0; /* kg/m^3 */ } // TODO: only for silicon (generalize)
      double get_phonon_energy() const { return 1.24 * (viennashe::physics::constants::q/20.0); /* Joule */ }
      double get_coupling_constant(viennashe::carrier_type_id) const { return 5e10 * viennashe::physics::constants::q; /* J/m */  }

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_   /* vecci: 0.48 */
                                                      : fit_p_;  /* vecci: 0.9 */
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };

    /** @brief Parameters for ionized impurity scattering using an isotropic fit, cf. papers by Jungemann */
    class ionized_impurity_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

      public:
        ionized_impurity_scattering_parameters() : base_type(false) {}

        double fit_factor(double NI, bool is_donor_dominated, viennashe::carrier_type_id ctype) const
        {
          if (ctype == viennashe::ELECTRON_TYPE_ID)
            return fit_factor_n(NI, is_donor_dominated);
          return fit_factor_p(NI, is_donor_dominated);
        }

      private:
        /** @brief Piecewise constant approximation of the fit factor according to Jungemann and Hong */
        double fit_factor_n(double NI, bool is_donor_dominated) const
        {
          double factor = 1.0;

          const long i = static_cast<long>(std::log(NI/100/100/100)/std::log(10.0));

          if (is_donor_dominated)
          {
           if (i <= 14) factor = 0.0479266894808;
           else if (i <= 15) factor = 0.0143814762161;
           else if (i <= 16) factor = 0.163014476723;
           else if (i <= 17) factor = 1.27621982607;
           else if (i <= 18) factor = 2.07841863841;
           else if (i <= 19) factor = 2.66848494897;
           else if (i <= 20) factor = 8.69553015708;
           else factor = 71.4713817122;
          }
          else
          {
           if (i <= 14) factor = 0.0135170431921;
           else if (i <= 15) factor = 0.110716650784;
           else if (i <= 16) factor = 0.499715004949;
           else if (i <= 17) factor = 1.4012436912;
           else if (i <= 18) factor = 1.24925936014;
           else if (i <= 19) factor = 1.0073534599;
           else if (i <= 20) factor = 2.84325491032;
           else factor = 22.5197987554;
          }

          return factor;
        }

        double fit_factor_p(double NI, bool is_donor_dominated) const
        {
          double factor = 1.0;

          const long i = static_cast<long>(std::log(NI/100/100/100)/std::log(10.0));


          if (is_donor_dominated)
          {

           if (i <= 14) factor = 0.0001;
           else if (i <= 15) factor = 0.0229086958854;
           else if (i <= 16) factor = 0.106259676639;
           else if (i <= 17) factor = 0.492336240353;
           else if (i <= 18) factor = 0.529652548279;
           else if (i <= 19) factor = 0.336197331939; // TODO: check this
           else if (i <= 20) factor = 0.512045950387;
           else factor = 3.01287648798;

          }
          else
          {

           if (i <= 14) factor =  1.1566304323;
           else if (i <= 15) factor = 0.324429154352; // TODO: check this
           else if (i <= 16) factor = 0.253416326007; // TODO: check this
           else if (i <= 17) factor = 0.718942981114; // TODO: check this
           else if (i <= 18) factor = 1.22458439328;
           else if (i <= 19) factor = 1.12207235068;
           else if (i <= 20) factor = 1.69835752347;
           else factor = 9.99694686735;

          }

          return factor;
        }

    };


    class impact_ionization_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      impact_ionization_scattering_parameters() : base_type(false), fit_n_(1.0), fit_p_(1.0) {}

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_ : fit_p_;
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };


    class trapped_charge_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      trapped_charge_scattering_parameters() : base_type(true), fit_n_(1.0), fit_p_(1.0)  {}

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_ : fit_p_;
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };

    class fixed_charge_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      fixed_charge_scattering_parameters() : base_type(false), fit_n_(1.0), fit_p_(1.0) {}

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_ : fit_p_;
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };


    class surface_scattering_parameters : public scattering_parameter_base
    {
      typedef scattering_parameter_base   base_type;

    public:
      surface_scattering_parameters() : base_type(false), fit_n_(1.0), fit_p_(1.0) {}

       // TODO: ... this is for electrons
      double first_factor(viennashe::carrier_type_id) const  { return 122238.00206847844; /* s V/m */ }

       // TODO: ... this is for electrons
      double second_factor(viennashe::carrier_type_id) const  { return 2.5541436127255013; /* s V^1/3  *  m^-1/3 */ }

       // TODO: ... this is for electrons
      double third_factor(viennashe::carrier_type_id) const  { return 1.485497231030198e-12; /* m^2 * V^-2 s^-1 */ }

      double cutoff_distance() const { return 10e-9 ; /* m */ }

      double get_fit_factor(viennashe::carrier_type_id ctype) const
      {
        return (ctype == viennashe::ELECTRON_TYPE_ID) ? fit_n_ : fit_p_;
      }

      void set_fit_factor(viennashe::carrier_type_id ctype, double value)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          fit_n_ = value;
        else
          fit_p_ = value;
      }

    private:
      double fit_n_;
      double fit_p_;
    };

#ifdef VIENNASHE_USE_DISABLED_CODE
    struct surface_roughness_scattering_parameters
    {

      // TODO: ... this is for electrons
      double get_rms_height(viennashe::carrier_type_id) const  { return 15e-10; /* m */ }

      // TODO: ... this is for electrons
      double get_correlation_length(viennashe::carrier_type_id) const  { return 29e-10; /* m */ }

      double get_fit_factor() const { return 1.0; }
    };


    struct surface_acoustic_phonon_scattering_parameters
    {
      double get_mass_density(viennashe::carrier_type_id) const { return 2330.0; /* kg/m^3 */ } // TODO: only for silicon (generalize)

      double get_longitudinal_sound_velocity(viennashe::carrier_type_id) const { return 9040.0; /* m/s */  }

      double get_deformation_potential(viennashe::carrier_type_id ctype) const
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return 8.9 * viennashe::physics::constants::q; /* Joule */
        else
          return 5.12 * viennashe::physics::constants::q;
      }

      // TODO: holes ... this is for electrons
      double get_potential_well_strength(viennashe::carrier_type_id) const  { return 13.1; }

      double get_fit_factor() const { return 1.0; }
    };
#endif


    //////////////////////////////////////


    /** @brief A configuration class for scattering mechanisms. Enable or disable scattering mechanisms here */
    class scatter_config
    {
      public:
        scatter_config() : scatter_ee_(false) {}

        // Scattering:
        /** @brief Returns the parameters for acoustic phonon scattering. Const-version. */
        acoustic_phonon_scattering_parameters const & acoustic_phonon() const { return scatter_ac_; }
        /** @brief Returns the parameters for acoustic phonon scattering. Non-const-version. */
        acoustic_phonon_scattering_parameters       & acoustic_phonon()       { return scatter_ac_; }

        /** @brief Returns the parameters for optical phonon scattering. Const-version. */
        optical_phonon_scattering_parameters const & optical_phonon() const { return scatter_op_; }
        /** @brief Returns the parameters for optical phonon scattering. Non-const-version. */
        optical_phonon_scattering_parameters       & optical_phonon()       { return scatter_op_; }

        /** @brief Returns the parameters for ionized impurity scattering. Const-version. */
        ionized_impurity_scattering_parameters const & ionized_impurity() const { return scatter_imp_; }
        /** @brief Returns the parameters for ionized impurity scattering. Non-const-version. */
        ionized_impurity_scattering_parameters       & ionized_impurity()       { return scatter_imp_; }

        /** @brief Returns the parameters for impact ionization scattering. Const-version. */
        impact_ionization_scattering_parameters const & impact_ionization() const { return scatter_impact_; }
        /** @brief Returns the parameters for impact ionization scattering. Non-const-version. */
        impact_ionization_scattering_parameters       & impact_ionization()       { return scatter_impact_; }

        /** @brief Returns true if electron-electron scattering is activated */
        bool electron_electron() const { return scatter_ee_; }
        /** @brief Enables/Disables electron-electron scattering */
        void electron_electron(bool b) { scatter_ee_ = b; }

        /** @brief Returns the parameters for fixed charge scattering. Const-version. */
        fixed_charge_scattering_parameters const & fixed_charge() const { return scatter_fixed_charge_; }
        /** @brief Returns the parameters for fixed charge scattering. Non-const-version. */
        fixed_charge_scattering_parameters       & fixed_charge()       { return scatter_fixed_charge_; }

        /** @brief Returns the parameters for fixed charge scattering. Const-version. */
        trapped_charge_scattering_parameters const & trapped_charge() const { return scatter_trapped_charge_; }
        /** @brief Returns the parameters for fixed charge scattering. Non-const-version. */
        trapped_charge_scattering_parameters       & trapped_charge()       { return scatter_trapped_charge_; }

        /** @brief Returns the parameters for surface scattering. Const-version. */
        surface_scattering_parameters const & surface() const { return scatter_surf_; }
        /** @brief Returns the parameters for surface scattering. Non-const-version. */
        surface_scattering_parameters       & surface()       { return scatter_surf_; }

      private:
        acoustic_phonon_scattering_parameters   scatter_ac_;
        optical_phonon_scattering_parameters    scatter_op_;
        ionized_impurity_scattering_parameters  scatter_imp_;
        impact_ionization_scattering_parameters scatter_impact_;
        bool scatter_ee_;     // electron-electron scattering, TODO: Add parameter class
        fixed_charge_scattering_parameters      scatter_fixed_charge_;
        trapped_charge_scattering_parameters    scatter_trapped_charge_;
        surface_scattering_parameters           scatter_surf_;
    };



  }


}

#endif
