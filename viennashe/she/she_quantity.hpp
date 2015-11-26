#ifndef VIENNASHE_SHE_SHE_QUANTITY_HPP
#define VIENNASHE_SHE_SHE_QUANTITY_HPP

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


/** @file viennashe/she/she_quantity.hpp
    @brief Defines a SHE quantity in (x, H)-space for use within the solvers of ViennaSHE
*/

#include <vector>
#include <string>

#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"

namespace viennashe
{
  namespace she
  {
    /** @brief Returns the number of even spherical harmonics up to (including) order L */
    inline long even_unknowns_on_node(long L)
    {
      long num = 0;
      for (long l=0; l<=L; l+=2)
        num += 2*l + 1;

      return num;
    }

    /** @brief Returns the number of odd spherical harmonics up to (including) order L */
    inline long odd_unknowns_on_node(long L)
    {
      long num = 0;
      for (long l=1; l<=L; l+=2)
        num += 2*l + 1;

      return num;
    }

    /** @brief General representation of any solver quantity defined on two different element types (e.g. vertices and edges) in an augmented (x, H) space */
    template<typename AssociatedT1,
             typename AssociatedT2,
             typename ValueT = double>
    class unknown_she_quantity
    {
      std::size_t get_id(AssociatedT1 const & elem) const { return std::size_t(viennagrid_index_from_element_id(elem)); }
      //std::size_t get_id(AssociatedT2 const & elem) const { return static_cast<std::size_t>(elem.id().get()); }

      public:
        typedef ValueT         value_type;

        unknown_she_quantity() {}  // to fulfill default constructible concept!

        unknown_she_quantity(std::string const & quan_name,
                             viennashe::carrier_type_id ctype,
                             equation_id quan_equation)
          : name_(quan_name),
            ctype_(ctype),
            equation_(quan_equation),
            log_damping_(false)
        {}

        void resize(std::size_t num_values_1, std::size_t num_values_2)
        {
          boundary_types1_.resize(num_values_1, BOUNDARY_NONE);
          boundary_types2_.resize(num_values_2, BOUNDARY_NONE);
          spatial_mask1_.resize(num_values_1, false);
          spatial_mask2_.resize(num_values_2, false);
          expansion_order_adaption_.resize(num_values_1);
          bandedge_shift1_.resize(num_values_1);
          bandedge_shift2_.resize(num_values_2);
        }

        void resize(std::size_t num_values_1, std::size_t num_values_2, std::size_t size_index_H)
        {
          resize(num_values_1, num_values_2);
          values1_.resize(num_values_1 * size_index_H);
          values2_.resize(num_values_2 * size_index_H);
          boundary_values1_.resize(num_values_1 * size_index_H);
          boundary_values2_.resize(num_values_2 * size_index_H);
          defined_but_unknown_mask1_.resize(num_values_1 * size_index_H, false);
          defined_but_unknown_mask2_.resize(num_values_2 * size_index_H, false);
          unknowns_indices1_.resize(num_values_1 * size_index_H, -1);
          unknowns_indices2_.resize(num_values_2 * size_index_H, -1);
          expansion_order1_.resize(num_values_1 * size_index_H);
          expansion_order2_.resize(num_values_2 * size_index_H);
          values_H_.resize(size_index_H);
        }

        std::string get_name() const { return name_; }

        equation_id get_equation() const { return equation_; }

        ValueT const * get_values(AssociatedT1 const & elem, std::size_t index_H) const { return &(values1_.at(array_index(get_id(elem), index_H)).at(0)); }
        //ValueT const * get_values(AssociatedT2 const & elem, std::size_t index_H) const { return &(values2_.at(array_index(get_id(elem), index_H)).at(0)); }

        void set_values(AssociatedT1 const & elem, std::size_t index_H, ValueT const * values)
        {
          for (std::size_t i=0; i < this->get_unknown_num(elem, index_H); ++i)
            values1_.at(array_index(get_id(elem), index_H)).at(i) = values[i];
        }
        /*void set_values(AssociatedT2 const & elem, std::size_t index_H, ValueT const * values)
        {
          for (std::size_t i=0; i < this->get_unknown_num(elem, index_H); ++i)
            values2_.at(array_index(get_id(elem), index_H)).at(i) = values[i];
        }*/

        // Dirichlet and Neumann
        ValueT get_boundary_value(AssociatedT1 const & elem, std::size_t index_H) const         { return boundary_values1_.at(array_index(get_id(elem), index_H));         }
        //ValueT get_boundary_value(AssociatedT2 const & elem, std::size_t index_H) const         { return boundary_values2_.at(array_index(get_id(elem), index_H));         }

        void   set_boundary_value(AssociatedT1 const & elem, std::size_t index_H, ValueT value) {        boundary_values1_.at(array_index(get_id(elem), index_H)) = value; }
        //void   set_boundary_value(AssociatedT2 const & elem, std::size_t index_H, ValueT value) {        boundary_values2_.at(array_index(get_id(elem), index_H)) = value; }

        boundary_type_id get_boundary_type(AssociatedT1 const & elem) const                   { return boundary_types1_.at(get_id(elem));         }
        //boundary_type_id get_boundary_type(AssociatedT2 const & elem) const                   { return boundary_types2_.at(get_id(elem));         }

        void             set_boundary_type(AssociatedT1 const & elem, boundary_type_id value) {        boundary_types1_.at(get_id(elem)) = value; }
        //void             set_boundary_type(AssociatedT2 const & elem, boundary_type_id value) {        boundary_types2_.at(get_id(elem)) = value; }

        // Unknown handling
        bool   get_unknown_mask(AssociatedT1 const & elem, std::size_t index_H) const       { return defined_but_unknown_mask1_.at(array_index(get_id(elem), index_H));         }
        //bool   get_unknown_mask(AssociatedT2 const & elem, std::size_t index_H) const       { return defined_but_unknown_mask2_.at(array_index(get_id(elem), index_H));         }

        void   set_unknown_mask(AssociatedT1 const & elem, std::size_t index_H, bool value) {        defined_but_unknown_mask1_.at(array_index(get_id(elem), index_H)) = value; }
        //void   set_unknown_mask(AssociatedT2 const & elem, std::size_t index_H, bool value) {        defined_but_unknown_mask2_.at(array_index(get_id(elem), index_H)) = value; }

        bool   get_unknown_mask(AssociatedT1 const & elem) const       { return spatial_mask1_.at(get_id(elem));         }
        //bool   get_unknown_mask(AssociatedT2 const & elem) const       { return spatial_mask2_.at(get_id(elem));         }

        void   set_unknown_mask(AssociatedT1 const & elem, bool value) {        spatial_mask1_.at(get_id(elem)) = value; }
        //void   set_unknown_mask(AssociatedT2 const & elem, bool value) {        spatial_mask2_.at(get_id(elem)) = value; }

        long   get_unknown_index(AssociatedT1 const & elem, std::size_t index_H) const       { return unknowns_indices1_.at(array_index(get_id(elem), index_H));         }
        //long   get_unknown_index(AssociatedT2 const & elem, std::size_t index_H) const       { return unknowns_indices2_.at(array_index(get_id(elem), index_H));         }

        void   set_unknown_index(AssociatedT1 const & elem, std::size_t index_H, long value) {        unknowns_indices1_.at(array_index(get_id(elem), index_H)) = value; }
        //void   set_unknown_index(AssociatedT2 const & elem, std::size_t index_H, long value) {        unknowns_indices2_.at(array_index(get_id(elem), index_H)) = value; }

        std::size_t  get_expansion_order(AssociatedT1 const & elem, std::size_t index_H) const        { return expansion_order1_.at(array_index(get_id(elem), index_H)); }
        //std::size_t  get_expansion_order(AssociatedT2 const & elem, std::size_t index_H) const        { return expansion_order2_.at(array_index(get_id(elem), index_H)); }

        void   set_expansion_order(AssociatedT1 const & elem, std::size_t index_H, std::size_t value)
        {
          expansion_order1_.at(array_index(get_id(elem), index_H)) = value;
          if (value > 0)
            values1_.at(array_index(get_id(elem), index_H)).resize(static_cast<std::size_t>(even_unknowns_on_node(static_cast<long>(this->get_expansion_order(elem, index_H)))));
          else
            values1_.at(array_index(get_id(elem), index_H)) = std::vector<ValueT>();
        }
        /*void   set_expansion_order(AssociatedT2 const & elem, std::size_t index_H, std::size_t value)
        {
          expansion_order2_.at(array_index(get_id(elem), index_H)) = value;
          if (value > 0)
            values2_.at(array_index(get_id(elem), index_H)).resize(static_cast<std::size_t>(odd_unknowns_on_node(static_cast<long>(this->get_expansion_order(elem, index_H)))));
          else
            values2_.at(array_index(get_id(elem), index_H)) = std::vector<ValueT>();
        }*/

        std::size_t  get_unknown_num(AssociatedT1 const & elem, std::size_t index_H) const
        {
          if ( this->get_unknown_index(elem, index_H) >= 0 )
            return static_cast<std::size_t>(even_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
          else
            return 0;
        }
        /*std::size_t  get_unknown_num(AssociatedT2 const & elem, std::size_t index_H) const
        {
          if ( this->get_unknown_index(elem, index_H) >= 0 )
            return  static_cast<std::size_t>(odd_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
          else
            return 0;
        }*/

        bool   get_expansion_adaption(AssociatedT1 const & elem) const   { return expansion_order_adaption_.at(get_id(elem)); }
        void   set_expansion_adaption(AssociatedT1 const & elem, bool b) {        expansion_order_adaption_.at(get_id(elem)) = b; }



        // total energy
        ValueT get_value_H(std::size_t index_H) const { return values_H_.at(index_H); }
        void   set_value_H(std::size_t index_H, ValueT value) { values_H_.at(index_H) = value; }

        std::size_t get_value_H_size() const { return values_H_.size(); }

        // band gap center
        ValueT get_bandedge_shift(AssociatedT1 const & elem) const { return bandedge_shift1_.at(get_id(elem)); }
        //ValueT get_bandedge_shift(AssociatedT2 const & elem) const { return bandedge_shift2_.at(get_id(elem)); }

        void set_bandedge_shift(AssociatedT1 const & elem, ValueT value) { bandedge_shift1_.at(get_id(elem)) = value; }
        //void set_bandedge_shift(AssociatedT2 const & elem, ValueT value) { bandedge_shift2_.at(get_id(elem)) = value; }

        ValueT get_kinetic_energy(AssociatedT1 const & elem, std::size_t index_H) const
        {
          return (ctype_ == ELECTRON_TYPE_ID) ? get_value_H(index_H) - bandedge_shift1_.at(get_id(elem))
                                              : bandedge_shift1_.at(get_id(elem)) - get_value_H(index_H);
        }
        /*ValueT get_kinetic_energy(AssociatedT2 const & elem, std::size_t index_H) const
        {
          return (ctype_ == ELECTRON_TYPE_ID) ? get_value_H(index_H) - bandedge_shift2_.at(get_id(elem))
                                              : bandedge_shift2_.at(get_id(elem)) - get_value_H(index_H);
        }*/

        carrier_type_id get_carrier_type_id() const { return ctype_; }

        /*
        std::size_t get_unknown_num() const
        {
          std::size_t num = 0;
          for (std::size_t i=0; i<unknowns_indices_.size(); ++i)
          {
            if (unknowns_indices_[i] >= 0)
              ++num;
          }
          return num;
        }*/

        // possible design flaws:
        bool get_logarithmic_damping() const { return log_damping_; }
        void set_logarithmic_damping(bool b) { log_damping_ = b; }

      private:
        std::size_t array_index(std::size_t element_id, std::size_t index_H) const { return element_id * values_H_.size() + index_H; }

        std::string                     name_;
        carrier_type_id                 ctype_;
        equation_id                     equation_;

        std::vector< std::vector<ValueT> >  values1_;
        std::vector< std::vector<ValueT> >  values2_;
        std::vector<std::size_t>            values1_offsets_;
        std::vector<std::size_t>            values2_offsets_;
        std::vector<boundary_type_id>       boundary_types1_;
        std::vector<boundary_type_id>       boundary_types2_;
        std::vector<ValueT>                 boundary_values1_;
        std::vector<ValueT>                 boundary_values2_;
        std::vector<bool>                   defined_but_unknown_mask1_;
        std::vector<bool>                   defined_but_unknown_mask2_;
        std::vector<bool>                   spatial_mask1_;
        std::vector<bool>                   spatial_mask2_;
        std::vector<long>                   unknowns_indices1_;
        std::vector<long>                   unknowns_indices2_;
        std::vector<std::size_t>            expansion_order1_;
        std::vector<std::size_t>            expansion_order2_;
        std::vector< bool >                 expansion_order_adaption_;

        energy_vector_type             values_H_;
        std::vector<ValueT>            bandedge_shift1_;
        std::vector<ValueT>            bandedge_shift2_;

        bool                           log_damping_;
    };


    /** @brief Common representation of any quantity associated with two objects of a certain type and an index */
    template<typename AssociatedT1, typename AssociatedT2,
             typename ValueT = double>
    class const_she_quantity
    {
      std::size_t get_id(AssociatedT1 const & elem) const { return static_cast<std::size_t>(elem.id().get()); }
      std::size_t get_id(AssociatedT2 const & elem) const { return static_cast<std::size_t>(elem.id().get()); }

      public:
        typedef ValueT         value_type;

        const_she_quantity(unknown_she_quantity<AssociatedT1, AssociatedT2> const & quan)
          : name_(quan.get_name()),
            values1_(quan.values1()),
            values2_(quan.values2()),
            expansion_order1_(quan.expansion_orders1()),
            expansion_order2_(quan.expansion_orders2())
        {}

        ValueT const * get_values(AssociatedT1 const & elem, std::size_t index_H) const { return &(values1_.at(get_id(elem)).at(index_H)); }
        ValueT const * get_values(AssociatedT2 const & elem, std::size_t index_H) const { return &(values2_.at(get_id(elem)).at(index_H)); }

        short  get_expansion_order(AssociatedT1 const & elem, std::size_t index_H) const        { return expansion_order1_.at(get_id(elem)).at(index_H);         }
        short  get_expansion_order(AssociatedT2 const & elem, std::size_t index_H) const        { return expansion_order2_.at(get_id(elem)).at(index_H);         }

        short  get_unknown_num(AssociatedT1 const & elem, std::size_t index_H) const { return even_unknowns_on_node(get_expansion_order(elem, index_H)); }
        short  get_unknown_num(AssociatedT2 const & elem, std::size_t index_H) const { return  odd_unknowns_on_node(get_expansion_order(elem, index_H)); }

        ValueT operator()(AssociatedT1 const & elem, std::size_t index_H) const { return values1_.at(get_id(elem)).at(index_H); }
        ValueT operator()(AssociatedT2 const & elem, std::size_t index_H) const { return values2_.at(get_id(elem)).at(index_H); }

        std::string name() const { return name_; }

      private:
        std::string                       name_;
        std::vector<std::vector<ValueT> > values1_;
        std::vector<std::vector<ValueT> > values2_;
        std::vector<std::vector<short> >  expansion_order1_;
        std::vector<std::vector<short> >  expansion_order2_;
    };

  } //namespace she
} //namespace viennashe

#endif
