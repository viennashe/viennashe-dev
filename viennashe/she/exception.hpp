#ifndef VIENNASHE_SHE_EXCEPTION_HPP
#define VIENNASHE_SHE_EXCEPTION_HPP

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
#include <stdexcept>
#include <sstream>

/** @file viennashe/she/exception.hpp
    @brief Provides the exceptions used inside the viennashe::she namespace.
*/

namespace viennashe
{
  namespace she
  {
    /** @brief Exception for the case that an invalid dispersion relation is specified */
    class unknown_dispersion_relation_exception : public std::runtime_error {
    public:
      unknown_dispersion_relation_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that no initial guess was/is specified */
    class no_init_guess_found_exception : public std::runtime_error {
    public:
      no_init_guess_found_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that a vertex has a coupling with itself */
    class coupled_vertices_equal_exception : public std::runtime_error {
    public:
      coupled_vertices_equal_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that any component encounters a divison by 0.0 */
    class division_by_zero : public std::runtime_error {
    public:
      division_by_zero(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that invalid expansion order is accessed
     *
     * Note that only even expansion orders are allowed on vertices, and odd expansion orders on edges.
     */
    class invalid_expansion_order_exception : public std::runtime_error
    {
      public:
        virtual const char* what() const throw() { return msg_.c_str(); }

        invalid_expansion_order_exception(std::string el,
                                          std::size_t l,
                                          long m) : std::runtime_error(el), element_string_(el), l_(l), m_(m)
        { this->fill_message(); }

        virtual ~invalid_expansion_order_exception() throw() {  }

      private:
        std::string element_string_;
        std::size_t l_;
        long m_;

        std::string msg_;

        void fill_message()
        {
          std::stringstream ss;
          ss << "* ViennaSHE: Invalid expansion order (" << l_ << ", " << m_ << ") on " << element_string_ << " accessed!";
          msg_ = ss.str();
        }
    };


    /** @brief Exception for the case that invalid expansion order is accessed
     *
     * Note that only even expansion orders are allowed on vertices, and odd expansion orders on edges.
     */
    class invalid_matrixelement_exception : public std::runtime_error
    {
      public:
        virtual const char* what() const throw() { return _msg.c_str(); }

        invalid_matrixelement_exception(std::string el, double value) : std::runtime_error(el), _value(value) { }

        virtual ~invalid_matrixelement_exception() throw() { }

      private:
        double _value;
        std::string _msg;

        void fill_message()
        {
          std::stringstream ss;
          ss << "* ViennaSHE: Invalid matrixelement '" << _value << "' found. ";
          _msg = ss.str();
        }
    };




    /** @brief Exception for the case that a macroscopic quantity is accessed, but the simulator has not yet been run. */
    class quantity_not_yet_available_exception : public std::runtime_error {
    public:
      quantity_not_yet_available_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception for the case that the total energy is smaller than the kinetic energy. */
    class total_energy_too_small_exception : public std::runtime_error {
    public:
      total_energy_too_small_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception for the case that neither electrons nor holes are selected for the simulation. */
    class no_carrier_type_id_specified_exception : public std::runtime_error {
    public:
      no_carrier_type_id_specified_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception for the case that neither electrons nor holes are selected for the simulation. */
    class she_simulator_does_not_accept_drift_diffusion_only_exception : public std::runtime_error {
    public:
      she_simulator_does_not_accept_drift_diffusion_only_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that traps are enabled without a bipolar SHE simulation. */
    class she_simulator_requires_bipolar_solution_for_traps : public std::runtime_error {
    public:
      she_simulator_requires_bipolar_solution_for_traps(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that neither electrons nor holes are selected for the simulation. */
    class negative_integration_interval_length_exception : public std::runtime_error {
    public:
      negative_integration_interval_length_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception for the case that a scattering term is Inf, NaN or causes an invalid system matrix entry. */
    class quantum_correction_with_newton_not_supported_exception : public std::runtime_error {
    public:
      quantum_correction_with_newton_not_supported_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception for the case that a scattering term is Inf, NaN or causes an invalid system matrix entry. */
    class assembly_exception : public std::runtime_error {
    public:
      assembly_exception(std::string const & str) : std::runtime_error(str) {}
    };


    /** @brief Exception for the case that a scattering term is Inf, NaN or causes an invalid system matrix entry. */
    class invalid_scattering_term_exception : public viennashe::she::assembly_exception {
    public:
      invalid_scattering_term_exception(std::string const & str) : viennashe::she::assembly_exception(str) {}
    };


    /** @brief Exception for the case that adaptive SHE cannot deal with the provided configuration (currently holes or ee-scattering) */
    class adaptive_she_not_available_for_this_configuration_exception : public std::runtime_error {
    public:
      adaptive_she_not_available_for_this_configuration_exception(std::string const & str) : std::runtime_error(str) {}
    };

    /** @brief Exception thrown in case an unkown or unsupported carrier type is found */
    class unkown_carrier_type_exception : public std::runtime_error {
    public:
      unkown_carrier_type_exception(std::string const & str) : std::runtime_error(str) {}
    };

  }
}

#endif
