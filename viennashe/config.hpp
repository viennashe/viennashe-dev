#ifndef VIENNASHE_SHE_CONFIG_HPP
#define VIENNASHE_SHE_CONFIG_HPP

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
#include <memory>

#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"

#include "viennashe/exception.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/solvers/config.hpp"
#include "viennashe/she/scattering/config.hpp"

#include "viennashe/simulator_quantity.hpp"

/** @file viennashe/config.hpp
    @brief The SHE configuration class is defined here.
*/

namespace viennashe
{
  namespace detail
  {
    struct density_gradient_config
    {
      density_gradient_config() : lambda_(0.2), dirichlet_(-0.11), bnd_type_id_(viennashe::BOUNDARY_DIRICHLET)
      {
        robin_coeffs_.alpha = 0.0;
        robin_coeffs_.beta = 0.0;
      }

      density_gradient_config(double lambda, double alpha, double beta)
       : lambda_(lambda), dirichlet_(0.0), bnd_type_id_(viennashe::BOUNDARY_ROBIN)
      {
        robin_coeffs_.alpha = alpha;
        robin_coeffs_.beta  = beta;
      }

      double lambda() const   { return lambda_; }
      void   lambda(double l) { lambda_ = l; }

      double dirichlet_boundary_value() const   { return dirichlet_; }
      void   dirichlet_boundary_value(double d) { dirichlet_ = d;    }

      viennashe::boundary_type_id boundary_type() const { return bnd_type_id_; }

      void boundary_type(viennashe::boundary_type_id newid) { bnd_type_id_ = newid; }

      viennashe::robin_boundary_coefficients<double>   robin_coeffs() const { return robin_coeffs_; }
      viennashe::robin_boundary_coefficients<double> & robin_coeffs()       { return robin_coeffs_; }

    private:
      double lambda_;
      viennashe::robin_boundary_coefficients<double> robin_coeffs_;
      double dirichlet_;
      viennashe::boundary_type_id bnd_type_id_;
    };

  } // namespace detail

  /** @brief Provides IDs for the dispersion relations */
  struct dispersion_relation_ids
  {
      enum { parabolic_dispersion,        // parabolic dispersion for holes
             modena_dispersion,           // modena dispersion for holes
             ext_vecchi_dispersion        // extended Vecchi model, Jin et al. TED 2011
           };
  };


  /** @brief The boundary condition configuration for SHE   */
  class she_boundary_conditions_config
  {
  public:
    she_boundary_conditions_config() : she_bnd_id_(viennashe::BOUNDARY_GENERATION_RECOMBINATION), tau_(1e-16) {}

    /** @brief Returns the type of boundary conditions used by SHE */
    viennashe::boundary_type_id type() const { return she_bnd_id_; }

    /** @brief Setter for the type of boundary conditions used by SHE */
    void type(viennashe::boundary_type_id new_id)
    {
      if (new_id == viennashe::BOUNDARY_DIRICHLET || new_id == viennashe::BOUNDARY_GENERATION_RECOMBINATION)
        she_bnd_id_ = new_id;
      else
        throw invalid_boundary_condition_exception("Encountered invalid boundary type for SHE. Only BOUNDARY_DIRICHLET and BOUNDARY_GENERATION_RECOMBINATION supported!");
    }

    /** @brief Returns the recombination rate used in a Robin boundary condition */
    double generation_recombination_rate() const { return tau_; }
    /** @brief Sets the recombination rate used in a Robin boundary condition
     * @param new_tau Must be positive and non-negative. Unit: seconds
     */
    void   generation_recombination_rate(double new_tau) { if (new_tau > 0) tau_ = new_tau; }

  private:
    viennashe::boundary_type_id she_bnd_id_;
    double tau_;  // generation/recombination rate (0.1 fs). Might need further calibration. The shorter the rate, the stronger is the distribution function clamped to the contact
  };



  /** @brief The main SHE configuration class. To be adjusted by the user for his/her needs. */
  class config
    : public dispersion_relation_ids
  {
    public:
      typedef viennashe::solvers::linear_solver_config      linear_solver_config_type;
      typedef viennashe::solvers::nonlinear_solver_config   nonlinear_solver_config_type;
      typedef viennashe::physics::dispersion_proxy          dispersion_relation_type;


      config() :
       with_electrons_( true ),
       electron_equation_id_(EQUATION_CONTINUITY),
       with_holes_( true ),
       hole_equation_id_(EQUATION_CONTINUITY),
       with_traps_( false ),
       with_trap_selfconsistency_( false ),
       dispersion_relation_electrons_( new viennashe::physics::modena_dispersion<viennashe::materials::si>(viennashe::ELECTRON_TYPE_ID) ),
       dispersion_relation_holes_( new viennashe::physics::modena_dispersion<viennashe::materials::si>(viennashe::HOLE_TYPE_ID) ),
       she_discretization_id_(SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF),
       she_scaling_id_(SHE_SCALING_KINETIC_ENERGY),
       L_max_(1), adaptive_expansions_(false),
       energy_spacing_(31.0 * viennashe::physics::constants::q / 1000.0),
       min_kinetic_energy_range_electrons_(viennashe::physics::constants::q),
       min_kinetic_energy_range_holes_(viennashe::physics::constants::q),
       max_kinetic_energy_range_electrons_(5.0*viennashe::physics::constants::q),  // ~5 eV is the work function for most semiconductors
       max_kinetic_energy_range_holes_(5.0*viennashe::physics::constants::q),      // ~5 eV is the work function for most semiconductors
       use_h_transformation_(true),
       use_quantum_correction_(false),
       scatter_conf_(),
       linear_solver_conf_(),
       nonlinear_solver_conf_(),
       with_hde_(false),
       with_quantum_correction_(false),
       time_step_size_(0),
       she_boundary_conf_(),
       dg_config_electrons_(0.2, -5.5e-4, -7.5e-5),
       dg_config_holes_(0.22, -4.2e-4, 8.3e-5)
       {}

      config(config const & other) :
       with_electrons_( other.with_electrons_ ),
       electron_equation_id_(other.electron_equation_id_),
       with_holes_( other.with_holes_ ),
       hole_equation_id_(other.hole_equation_id_),
       with_traps_( other.with_traps_ ),
       with_trap_selfconsistency_(other.with_trap_selfconsistency_),
       dispersion_relation_electrons_(other.dispersion_relation_electrons_->clone()) ,
       dispersion_relation_holes_(other.dispersion_relation_holes_->clone()) ,
       she_discretization_id_(other.she_discretization_id_),
       she_scaling_id_(other.she_scaling_id_),
       L_max_(other.L_max_), adaptive_expansions_(other.adaptive_expansions_),
       energy_spacing_(other.energy_spacing_),
       min_kinetic_energy_range_electrons_(other.min_kinetic_energy_range_electrons_),
       min_kinetic_energy_range_holes_(other.min_kinetic_energy_range_holes_),
       max_kinetic_energy_range_electrons_(other.max_kinetic_energy_range_electrons_),
       max_kinetic_energy_range_holes_(other.max_kinetic_energy_range_holes_),
       use_h_transformation_(other.use_h_transformation_),
       use_quantum_correction_(other.use_quantum_correction_),
       scatter_conf_(other.scatter_conf_),
       linear_solver_conf_(other.linear_solver_conf_),
       nonlinear_solver_conf_(other.nonlinear_solver_conf_),
       with_hde_(other.with_hde_),
       with_quantum_correction_(other.with_quantum_correction_),
       time_step_size_(other.time_step_size_),
       she_boundary_conf_(other.she_boundary_conf_),
       dg_config_electrons_(other.dg_config_electrons_),
       dg_config_holes_(other.dg_config_holes_)
       {}

      void operator=(config const & other)
      {
        with_electrons_ = other.with_electrons_;
        with_holes_ = other.with_holes_;
        with_traps_ = other.with_traps_;
        with_trap_selfconsistency_ = other.with_trap_selfconsistency_;
//        dispersion_relation_electrons_ = std::auto_ptr< viennashe::physics::dispersion_base >(other.dispersion_relation_electrons_->clone());
//        dispersion_relation_holes_ = std::auto_ptr< viennashe::physics::dispersion_base >(other.dispersion_relation_holes_->clone());
        dispersion_relation_electrons_ = std::unique_ptr< viennashe::physics::dispersion_base >(other.dispersion_relation_electrons_->clone());
        dispersion_relation_holes_ = std::unique_ptr< viennashe::physics::dispersion_base >(other.dispersion_relation_holes_->clone());
        she_discretization_id_ = other.she_discretization_id_;
        she_scaling_id_        = other.she_scaling_id_;
        L_max_ = other.L_max_;
        adaptive_expansions_ = other.adaptive_expansions_;
        energy_spacing_ = other.energy_spacing_;
        min_kinetic_energy_range_electrons_ = other.min_kinetic_energy_range_electrons_;
        min_kinetic_energy_range_holes_ = other.min_kinetic_energy_range_holes_;
        max_kinetic_energy_range_electrons_ = other.max_kinetic_energy_range_electrons_;
        max_kinetic_energy_range_holes_ = other.max_kinetic_energy_range_holes_;
        use_h_transformation_ = other.use_h_transformation_;
        use_quantum_correction_ = other.use_quantum_correction_;
        scatter_conf_ = other.scatter_conf_;
        linear_solver_conf_ = other.linear_solver_conf_;
        nonlinear_solver_conf_ = other.nonlinear_solver_conf_;
        with_hde_ = other.with_hde_;
        with_quantum_correction_ = other.with_quantum_correction_;
        time_step_size_ = other.time_step_size_;
        she_boundary_conf_ = other.she_boundary_conf_;
        dg_config_electrons_ = other.dg_config_electrons_;
        dg_config_holes_ = other.dg_config_holes_;
      }

      /////// Polarity //////////

      /** @brief Returns true if electrons are considered in the simulation */
      bool with_electrons() const { return with_electrons_; }
      /** @brief Activates are deactivates electrons in the simulation */
      void with_electrons(bool b) { with_electrons_ = b; }

      equation_id get_electron_equation() const { return electron_equation_id_; }
      void        set_electron_equation(equation_id equ_id) { electron_equation_id_ = equ_id; }

      /** @brief Returns true if holes are considered in the simulation */
      bool with_holes() const { return with_holes_; }
      /** @brief Activates are deactivates holes in the simulation */
      void with_holes(bool b) { with_holes_ = b; }

      equation_id get_hole_equation() const { return hole_equation_id_; }
      void        set_hole_equation(equation_id equ_id) { hole_equation_id_ = equ_id; }


      /** @brief Returns true if traps are considered in the simulation */
      bool with_traps() const { return with_traps_; }
      /** @brief Activates or deactivates traps in the simulation */
      void with_traps(bool b) { with_traps_ = b; }
      /** @brief Returns true if traps are considered self-consistently in the simulation */
      bool with_trap_selfconsistency() const { return with_trap_selfconsistency_; }
      /** @brief Activates or deactivates trap self-consistency in the simulation */
      void with_trap_selfconsistency(bool b) { with_trap_selfconsistency_ = b; }



      /** @brief Returns the dispersion relation for electrons */
      viennashe::physics::dispersion_proxy dispersion_relation_electrons() const
      {
        return viennashe::physics::dispersion_proxy(dispersion_relation_electrons_.get());
      }

      /** @brief Returns the dispersion relation for holes */
      viennashe::physics::dispersion_proxy dispersion_relation_holes() const
      {
        return viennashe::physics::dispersion_proxy(dispersion_relation_holes_.get());
      }

      /** @brief Returns the dispersion relation for electrons */
      viennashe::physics::dispersion_proxy dispersion_relation(viennashe::carrier_type_id ctype) const
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return viennashe::physics::dispersion_proxy(dispersion_relation_electrons_.get());
        else if (ctype == viennashe::HOLE_TYPE_ID)
          return viennashe::physics::dispersion_proxy(dispersion_relation_holes_.get());
        else
          throw viennashe::carrier_type_not_supported_exception("In dispersion_relation");
      }

      /** @brief Sets a new dispersion relation for electrons
       *
       * @param dispersion_id     Identifier of the dispersion relation, @see dispersion_relation_ids
       */
      void dispersion_relation_electrons(long dispersion_id)
      {
        viennashe::physics::dispersion_base * disp_ptr
             = dispersion_relation_impl<viennashe::materials::si>(viennashe::ELECTRON_TYPE_ID, dispersion_id);

        if (disp_ptr == NULL)
        {
          std::stringstream ss;
          ss << "No electron dispersion found with id: " << dispersion_id << std::endl;
          throw she::unknown_dispersion_relation_exception(ss.str());
        }
        dispersion_relation_electrons_ = std::unique_ptr< viennashe::physics::dispersion_base >( disp_ptr );

    //    dispersion_relation_electrons_ = std::auto_ptr< viennashe::physics::dispersion_base >( disp_ptr );
      }

      /** @brief Sets a new dispersion relation for electrons
       *
       * @param dispersion_id     Identifier of the dispersion relation, @see dispersion_relation_ids
       * @param ctype             Identifier of the carrier type
       */
      void dispersion_relation(long dispersion_id, viennashe::carrier_type_id ctype)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          dispersion_relation_electrons(dispersion_id);
        else
          dispersion_relation_holes(dispersion_id);
      }

      /** @brief Sets a new dispersion relation for electrons and holes
       *
       * @param name     String identifier of the dispersion relation
       */
      void dispersion_relation(const std::string name)
      {
        this->dispersion_relation_electrons(this->get_dispersion_relation_id_impl(name));
        this->dispersion_relation_holes(this->get_dispersion_relation_id_impl(name));
      }

      /** @brief Sets a new dispersion relation for electrons
       *
       * @param name     String identifier of the dispersion relation
       */
      void dispersion_relation_electrons(const std::string name)
      {
        this->dispersion_relation_electrons(this->get_dispersion_relation_id_impl(name));
      }


      /** @brief Sets a new dispersion relation for holes
       *
       * @param dispersion_id     Identifier of the dispersion relation, @see dispersion_relation_ids
       */
      void dispersion_relation_holes(long dispersion_id)  //TODO: Think about incorporation of materials here
      {
        viennashe::physics::dispersion_base * disp_ptr
             = dispersion_relation_impl<viennashe::materials::si>(viennashe::HOLE_TYPE_ID, dispersion_id);

        if (disp_ptr == NULL)
        {
          std::stringstream ss;
          ss << "No hole dispersion found with id: " << dispersion_id << std::endl;
          throw she::unknown_dispersion_relation_exception(ss.str());
        }

        //dispersion_relation_holes_ = std::auto_ptr< viennashe::physics::dispersion_base >( disp_ptr );
        dispersion_relation_holes_ = std::unique_ptr< viennashe::physics::dispersion_base >( disp_ptr );

      }

      /** @brief Sets a new dispersion relation for holes
       *
       * @param name     String identifier of the dispersion relation
       */
      void dispersion_relation_holes(const std::string name)
      {
        this->dispersion_relation_holes(this->get_dispersion_relation_id_impl(name));
      }


      she_discretization_type_id she_discretization_type() const { return she_discretization_id_; }
      void she_discretization_type(she_discretization_type_id discretization_id) { she_discretization_id_ = discretization_id; }

      she_scaling_type_id she_scaling_type() const { return she_scaling_id_; }
      void she_scaling_type(she_scaling_type_id scaling_id) { she_scaling_id_ = scaling_id; }


      // Expansion order:
      /** @brief Returns the current maximum expansion order */
      long max_expansion_order() const { return L_max_; }

      /** @brief Sets a new maximum expansion order. For uniform expansions, new_L will be used all over the device. For adaptive SHE, adaption is stopped at new_L */
      void max_expansion_order(long new_L)
      {
        assert(new_L > 0);
        L_max_ = new_L;
      }

      /** @brief Returns the flag for the use of adaptive expansions */
      bool adaptive_expansions() const { return adaptive_expansions_; }
      /** @brief Sets the use of adaptive expansions */
      void adaptive_expansions(bool b) { adaptive_expansions_ = b; }

      //
      // Minimum kinetic energy range:
      //

      /** @brief Returns the minimum kinetic energy range for the selected carrier */
      double min_kinetic_energy_range(viennashe::carrier_type_id ctype) const
      {
        if (ctype == ELECTRON_TYPE_ID)
          return min_kinetic_energy_range_electrons_;
        else
          return min_kinetic_energy_range_holes_;
      }

      /** @brief Sets the minimum kinetic energy range for the selected carrier */
      void min_kinetic_energy_range(double e_new, viennashe::carrier_type_id ctype)
      {
        if (e_new > 0)
        {
          if (ctype == ELECTRON_TYPE_ID)
            min_kinetic_energy_range_electrons_ = e_new;
          else
            min_kinetic_energy_range_holes_ = e_new;
        }
      }


      // set both
      /** @brief Sets the minimum kinetic energy range in the conduction band for electrons and in the valence band for electrons */
      void min_kinetic_energy_range(double e_new)
      {
        if (e_new > 0)
        {
          min_kinetic_energy_range_electrons_ = e_new;
          min_kinetic_energy_range_holes_ = e_new;
        }
      }

      //
      // Maximum kinetic energy range:
      //

      /** @brief Returns the minimum kinetic energy range for the selected carrier */
      double max_kinetic_energy_range(viennashe::carrier_type_id ctype) const
      {
        if (ctype == ELECTRON_TYPE_ID)
          return max_kinetic_energy_range_electrons_;
        else
          return max_kinetic_energy_range_holes_;
      }

      /** @brief Sets the minimum kinetic energy range for the selected carrier */
      void max_kinetic_energy_range(double e_new, viennashe::carrier_type_id ctype)
      {
        if (e_new > 0)
        {
          if (ctype == ELECTRON_TYPE_ID)
            max_kinetic_energy_range_electrons_ = e_new;
          else
            max_kinetic_energy_range_holes_ = e_new;
        }
      }


      // set both
      /** @brief Sets the minimum kinetic energy range in the conduction band for electrons and in the valence band for electrons */
      void max_kinetic_energy_range(double e_new)
      {
        if (e_new > 0)
        {
          max_kinetic_energy_range_electrons_ = e_new;
          max_kinetic_energy_range_holes_ = e_new;
        }
      }



      /** @brief Returns the uniform energy spacing of discrete energies */
      double energy_spacing() const { return energy_spacing_; }
      /** @brief Sets a new discrete energy spacing. */
      void energy_spacing(double new_spacing) { energy_spacing_ = new_spacing; }

      /** @brief Returns whether the H-transformation is used */
      bool use_h_transformation() const { return use_h_transformation_; }
      /** @brief Sets whether the H-transformation should be used for the energy discretization */
      void use_h_transformation(bool b) { use_h_transformation_ = b; }


      //
      // scattering
      //

      /** @brief Returns the configuration object for scattering */
      viennashe::she::scatter_config       & scattering()       { return scatter_conf_; }
      /** @brief Returns the configuration object for scattering */
      viennashe::she::scatter_config const & scattering() const { return scatter_conf_; }


      //
      // linear solver
      //

      /** @brief Returns the configuration object for the linear solver */
      linear_solver_config_type       & linear_solver()       { return linear_solver_conf_; }
      /** @brief Returns the configuration object for the linear solver */
      linear_solver_config_type const & linear_solver() const { return linear_solver_conf_; }


      //
      // nonlinear solver
      //

      /** @brief Returns the configuration object for the nonlinear solver */
      nonlinear_solver_config_type       & nonlinear_solver()       { return nonlinear_solver_conf_; }
      /** @brief Returns the configuration object for the nonlinear solver */
      nonlinear_solver_config_type const & nonlinear_solver() const { return nonlinear_solver_conf_; }

      //
      // HDE
      //
      bool with_hde() const { return with_hde_; }
      void with_hde(bool v) { with_hde_ = v;    }

      //
      // Quantum Correction
      //

      bool quantum_correction() const { return use_quantum_correction_; }
      void quantum_correction(bool use_quantum_correction)  { use_quantum_correction_ = use_quantum_correction; }

      bool with_quantum_correction() const { return with_quantum_correction_; }
      void with_quantum_correction(bool b) { with_quantum_correction_ = b; }

      //
      // Density Gradient
      //

      detail::density_gradient_config const & density_gradient(viennashe::carrier_type_id ctype) const
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return dg_config_electrons_;
        else
          return dg_config_holes_;
      }

      detail::density_gradient_config & density_gradient(viennashe::carrier_type_id ctype)
      {
        if (ctype == viennashe::ELECTRON_TYPE_ID)
          return dg_config_electrons_;
        else
          return dg_config_holes_;
      }

      ////////////////
      double time_step_size() const   { return time_step_size_; }
      //void   time_step_size(double s) { assert(s >= 0 && bool("Time step size must not be negative!")); time_step_size_ = s; }

      ////////////////
      she_boundary_conditions_config const & she_boundary_conditions() const { return she_boundary_conf_; }
      she_boundary_conditions_config       & she_boundary_conditions()       { return she_boundary_conf_; }

      //
      // Insulator Distance
      //
      bool setup_insulator_distances()       const { return this->scattering().surface().enabled(); }

    private:

      template <typename MaterialTag>
      viennashe::physics::dispersion_base * dispersion_relation_impl(viennashe::carrier_type_id ctype, long dispersion_id)
      {
        switch(dispersion_id)
        {
          case parabolic_dispersion:
            return new viennashe::physics::parabolic_dispersion<MaterialTag>(ctype);
          case modena_dispersion:
            return new viennashe::physics::modena_dispersion<MaterialTag>(ctype);
          case ext_vecchi_dispersion:
            return new viennashe::physics::ext_vecchi_dispersion<MaterialTag>(ctype);
          default:
            return NULL;
        };
        return NULL;
      }

      long get_dispersion_relation_id_impl(std::string name)
      {
        // to upper to ease the comparison and to be case insensitive
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);

        if (name == "parabolic_dispersion")
        {
          return parabolic_dispersion;
        }
        else if (name == "modena_dispersion")
        {
          return modena_dispersion;
        }
        else if (name == "ext_vecchi_dispersion" || name == "fullband")
        {
          return ext_vecchi_dispersion;
        }
        else
        {
          return -1;
        }
      }

      bool with_electrons_;
      equation_id electron_equation_id_;
      bool with_holes_;
      equation_id hole_equation_id_;
      bool with_traps_;
      bool with_trap_selfconsistency_;
      std::unique_ptr< viennashe::physics::dispersion_base >  dispersion_relation_electrons_;
      std::unique_ptr< viennashe::physics::dispersion_base >  dispersion_relation_holes_;
      she_discretization_type_id  she_discretization_id_;
      she_scaling_type_id         she_scaling_id_;

      long L_max_;
      bool adaptive_expansions_;

      double energy_spacing_;
      double min_kinetic_energy_range_electrons_;
      double min_kinetic_energy_range_holes_;
      double max_kinetic_energy_range_electrons_;
      double max_kinetic_energy_range_holes_;
      bool   use_h_transformation_;

      bool use_quantum_correction_;
      viennashe::she::scatter_config         scatter_conf_;
      linear_solver_config_type        linear_solver_conf_;
      nonlinear_solver_config_type  nonlinear_solver_conf_;

      bool with_hde_;
      bool with_quantum_correction_;

      double time_step_size_;
      she_boundary_conditions_config she_boundary_conf_;

      detail::density_gradient_config dg_config_electrons_;
      detail::density_gradient_config dg_config_holes_;


  };

} // namespace viennashe

#endif
