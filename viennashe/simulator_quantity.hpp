#ifndef VIENNASHE_SIMULATOR_QUANTITY_HPP
#define VIENNASHE_SIMULATOR_QUANTITY_HPP

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


/** @file viennashe/simulator_quantity.hpp
    @brief Defines a generic simulator quantity for use within the solvers of ViennaSHE
*/

#include <vector>
#include <string>

#include "viennagrid/mesh/mesh.hpp"

#include "viennashe/forwards.h"

namespace viennashe
{

  /** @brief Holds a function per quantity, which returns the name of the respective quantity as std::string.  */
  namespace quantity
  {
    inline std::string potential()                               { return "Electrostatic potential"; }
    inline std::string electron_density()                        { return "Electron density"; }
    inline std::string hole_density()                            { return "Hole density"; }
    inline std::string density_gradient_electron_correction()    { return "DG electron correction potential"; }
    inline std::string density_gradient_hole_correction()        { return "DG hole correction potential"; }
    inline std::string lattice_temperature()                     { return "Lattice temperature"; }
    inline std::string electron_distribution_function()          { return "Electron distribution function"; }
    inline std::string hole_distribution_function()              { return "Hole distribution function"; }
  }

  /** @brief Common representation of any quantity associated with objects of a certain type.
   *
   * This is the minimum requirement one can have: Return value for a vertex/cell.
   */
  template<typename AssociatedT, typename ValueT = double>
  class const_quantity
  {
    public:
      typedef const_quantity<AssociatedT, ValueT> self_type;

      typedef ValueT          value_type;
      typedef AssociatedT     associated_type;
      typedef associated_type access_type;

      const_quantity(std::string quan_name,
                     std::size_t num_values,
                     value_type default_value = value_type())
        : name_(quan_name),
          values_(num_values, default_value) {}

      const_quantity(std::string quan_name,
                     std::vector<ValueT> const & values_array)
        : name_(quan_name),
          values_(values_array) {}

      const_quantity(self_type const & o) : name_(o.name_), values_(o.values_) { }
      void operator=(self_type const & o) { name_(o.name_); values_(o.values_); }

      ValueT get_value (associated_type const & elem) const { return values_.at(static_cast<std::size_t>(elem.id().get())); }
      ValueT at        (associated_type const & elem) const { return this->get_value(elem); }
      ValueT operator()(associated_type const & elem) const { return this->get_value(elem); }

      std::string name() const { return name_; }

    private:
      std::string         name_;
      std::vector<ValueT> values_;
  };

  template<typename AssociatedT, typename ValueT = double>
  class non_const_quantity
  {
    public:
      typedef non_const_quantity<AssociatedT, ValueT> self_type;

      typedef ValueT         value_type;
      typedef AssociatedT    associated_type;

      non_const_quantity(std::string quan_name,
                         std::size_t num_values,
                         value_type default_value = value_type())
        : name_(quan_name), values_(num_values, default_value) {}

      non_const_quantity(std::string quan_name,
                         std::vector<ValueT> const & values_array)
        : name_(quan_name), values_(values_array) {}

      non_const_quantity(self_type const & o) : name_(o.name_), values_(o.values_) { }
      void operator=(self_type const & o) { name_(o.name_); values_(o.values_); }


      ValueT get_value(associated_type const & elem)  const         { return values_.at(static_cast<std::size_t>(elem.id().get()));   }
      ValueT at       (associated_type const & elem)  const         { return this->get_value(elem);   }
      ValueT operator()(associated_type const & elem) const         { return this->get_value(elem);   }
      void   set_value(associated_type const & elem, ValueT value)  { values_.at(static_cast<std::size_t>(elem.id().get())) = value;  }

      std::string name() const { return name_; }

    private:
      std::string         name_;
      std::vector<ValueT> values_;
  };

  template <typename ValueT = double>
  struct robin_boundary_coefficients
  {
    ValueT alpha;
    ValueT beta;
  };

  template<typename AssociatedT, typename ValueT = double>
  class unknown_quantity
  {
    std::size_t get_id(AssociatedT const & elem) const { return static_cast<std::size_t>(elem.id().get()); }

    public:
      typedef unknown_quantity<AssociatedT, ValueT> self_type;

      typedef ValueT         value_type;
      typedef AssociatedT    associated_type;
      typedef robin_boundary_coefficients<ValueT>   robin_boundary_type;

      unknown_quantity() {}  // to fulfill default constructible concept!

      unknown_quantity(std::string const & quan_name,
                       equation_id quan_equation,
                       std::size_t num_values,
                       value_type default_value = value_type())
        : name_(quan_name),
          equation_(quan_equation),
          values_                  (num_values, default_value),
          boundary_types_          (num_values, BOUNDARY_NONE),
          boundary_values_         (num_values, default_value),
          boundary_values_alpha_   (num_values, 0),
          defined_but_unknown_mask_(num_values, false),
          unknowns_indices_        (num_values, -1),
          log_damping_(false)
      {}

      unknown_quantity(self_type const & o)
        : name_(o.name_),
          equation_(o.equation_),
          values_                  (o.values_),
          boundary_types_          (o.boundary_types_),
          boundary_values_         (o.boundary_values_),
          boundary_values_alpha_   (o.boundary_values_alpha_),
          defined_but_unknown_mask_(o.defined_but_unknown_mask_),
          unknowns_indices_        (o.unknowns_indices_),
          log_damping_             (o.log_damping_)
      { }

      void operator=(self_type const & o)
      {
          name_ = o.name_;
          equation_ = o.equation_;
          values_                  = o.values_;
          boundary_types_          = o.boundary_types_;
          boundary_values_         = o.boundary_values_;
          boundary_values_alpha_   = o.boundary_values_alpha_;
          defined_but_unknown_mask_= o.defined_but_unknown_mask_;
          unknowns_indices_        = o.unknowns_indices_;
          log_damping_             = o.log_damping_;
      }

      std::string get_name() const { return name_; }

      equation_id get_equation() const { return equation_; }

      ValueT get_value(associated_type const & elem) const         { return values_.at(get_id(elem));  }
      ValueT at       (associated_type const & elem) const         { return this->get_value(elem);        }
      ValueT operator()(associated_type const & elem) const        { return this->get_value(elem);  }

      void   set_value(associated_type const & elem, ValueT value) {        values_.at(get_id(elem)) = value; }
      void   set_value(associated_type const & elem, robin_boundary_type value) { (void)elem; (void)value; } //TODO: Fix this rather ugly hack which stems from set_boundary_for_material()

      // Dirichlet and Neumann
      ValueT get_boundary_value(associated_type const & elem) const         { return boundary_values_.at(get_id(elem));         }
      void   set_boundary_value(associated_type const & elem, ValueT value) {        boundary_values_.at(get_id(elem)) = value; }

      // Robin boundary conditions
      /** @brief Returns the coefficients alpha and beta for the Robin boundary condition   du/dn = beta - alpha u
        *
        *
        * @return    A pair holding (beta, alpha). .first returns beta, .second returns alpha
        */
      robin_boundary_type get_boundary_values(associated_type const & elem) const
      {
        robin_boundary_type ret;
        ret.alpha = boundary_values_alpha_.at(get_id(elem));
        ret.beta  = boundary_values_.at(get_id(elem));
        return ret;
      }

      /** @brief Setter function for Robin boundary conditions.
        *
        * Note that this is an overload rather than a new function 'set_boundary_values' in order to have a uniform setter interface.
        */
      void set_boundary_value(associated_type const & elem, robin_boundary_type const & values)
      {
        boundary_values_.at(get_id(elem))       = values.beta;
        boundary_values_alpha_.at(get_id(elem)) = values.alpha;
      }

      boundary_type_id get_boundary_type(associated_type const & elem) const                   { return boundary_types_.at(get_id(elem));         }
      void             set_boundary_type(associated_type const & elem, boundary_type_id value) {        boundary_types_.at(get_id(elem)) = value; }

      bool   get_unknown_mask(associated_type const & elem) const       { return defined_but_unknown_mask_.at(get_id(elem));         }
      void   set_unknown_mask(associated_type const & elem, bool value) {        defined_but_unknown_mask_.at(get_id(elem)) = value; }

      long   get_unknown_index(associated_type const & elem) const       { return unknowns_indices_.at(get_id(elem));         }
      void   set_unknown_index(associated_type const & elem, long value) {        unknowns_indices_.at(get_id(elem)) = value; }

      std::size_t get_unknown_num() const
      {
        std::size_t num = 0;
        for (std::size_t i=0; i<unknowns_indices_.size(); ++i)
        {
          if (unknowns_indices_[i] >= 0)
            ++num;
        }
        return num;
      }

      // possible design flaws:
      std::vector<ValueT>           const & values()                   const { return values_; }
      std::vector<boundary_type_id> const & boundary_types()           const { return boundary_types_; }
      std::vector<ValueT>           const & boundary_values()          const { return boundary_values_; }
      std::vector<bool>             const & defined_but_unknown_mask() const { return defined_but_unknown_mask_; }
      std::vector<long>             const & unknowns_indices()         const { return unknowns_indices_; }


      bool get_logarithmic_damping() const { return log_damping_; }
      void set_logarithmic_damping(bool b) { log_damping_ = b; }

    private:
      std::string                     name_;
      equation_id                    equation_;
      std::vector<ValueT>            values_;
      std::vector<boundary_type_id>  boundary_types_;
      std::vector<ValueT>            boundary_values_;
      std::vector<ValueT>            boundary_values_alpha_;
      std::vector<bool>              defined_but_unknown_mask_;
      std::vector<long>              unknowns_indices_;
      bool                           log_damping_;
  };


} //namespace viennashe

#endif
