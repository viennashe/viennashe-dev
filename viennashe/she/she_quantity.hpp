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
    template<typename ValueT = double>
    class unknown_she_quantity
    {
      std::size_t get_id(viennagrid_element_id elem) const { return std::size_t(viennagrid_index_from_element_id(elem)); }

      public:
        typedef ValueT         value_type;

        unknown_she_quantity() {}  // to fulfill default constructible concept!

        unknown_she_quantity(viennagrid_mesh spatial_mesh,
                             viennagrid_dimension density_dimension,
                             viennagrid_dimension flux_dimension,
                             std::string const & quan_name,
                             viennashe::carrier_type_id ctype,
                             equation_id quan_equation)
          : mesh_(spatial_mesh),
            density_dim_(density_dimension),
            flux_dim_(flux_dimension),
            name_(quan_name),
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

        ValueT const * get_values(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
            return &(values1_.at(array_index(get_id(elem), index_H)).at(0));
          else if (topo_dim == flux_dim_)
            return &(values2_.at(array_index(get_id(elem), index_H)).at(0));
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_values()");
        }

        void set_values(viennagrid_element_id elem, std::size_t index_H, ValueT const * values)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            for (std::size_t i=0; i < this->get_unknown_num(elem, index_H); ++i)
              values1_.at(array_index(get_id(elem), index_H)).at(i) = values[i];
          }
          else if (topo_dim == flux_dim_)
          {
            for (std::size_t i=0; i < this->get_unknown_num(elem, index_H); ++i)
              values2_.at(array_index(get_id(elem), index_H)).at(i) = values[i];
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::set_values()");
        }


        // Dirichlet and Neumann
        ValueT get_boundary_value(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return boundary_values1_.at(array_index(get_id(elem), index_H));
          }
          else if (topo_dim == flux_dim_)
          {
            return boundary_values2_.at(array_index(get_id(elem), index_H));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::set_values()");
        }

        void set_boundary_value(viennagrid_element_id elem, std::size_t index_H, ValueT value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            boundary_values1_.at(array_index(get_id(elem), index_H)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            boundary_values2_.at(array_index(get_id(elem), index_H)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::set_values()");
        }

        boundary_type_id get_boundary_type(viennagrid_element_id elem) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return boundary_types1_.at(get_id(elem));
          }
          else if (topo_dim == flux_dim_)
          {
            return boundary_types2_.at(get_id(elem));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_boundary_type(viennagrid_element_id elem, boundary_type_id value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            boundary_types1_.at(get_id(elem)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            boundary_types2_.at(get_id(elem)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        // Unknown handling
        bool get_unknown_mask(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return defined_but_unknown_mask1_.at(array_index(get_id(elem), index_H));
          }
          else if (topo_dim == flux_dim_)
          {
            return defined_but_unknown_mask2_.at(array_index(get_id(elem), index_H));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_unknown_mask(viennagrid_element_id elem, std::size_t index_H, bool value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            defined_but_unknown_mask1_.at(array_index(get_id(elem), index_H)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            defined_but_unknown_mask2_.at(array_index(get_id(elem), index_H)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        bool get_unknown_mask(viennagrid_element_id elem) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return spatial_mask1_.at(get_id(elem));
          }
          else if (topo_dim == flux_dim_)
          {
            return spatial_mask2_.at(get_id(elem));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_unknown_mask(viennagrid_element_id elem, bool value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            spatial_mask1_.at(get_id(elem)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            spatial_mask2_.at(get_id(elem)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        long get_unknown_index(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return unknowns_indices1_.at(array_index(get_id(elem), index_H));
          }
          else if (topo_dim == flux_dim_)
          {
            return unknowns_indices2_.at(array_index(get_id(elem), index_H));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_unknown_index(viennagrid_element_id elem, std::size_t index_H, long value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            unknowns_indices1_.at(array_index(get_id(elem), index_H)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            unknowns_indices2_.at(array_index(get_id(elem), index_H)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        std::size_t get_expansion_order(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return expansion_order1_.at(array_index(get_id(elem), index_H));
          }
          else if (topo_dim == flux_dim_)
          {
            return expansion_order2_.at(array_index(get_id(elem), index_H));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_expansion_order(viennagrid_element_id elem, std::size_t index_H, std::size_t value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            expansion_order1_.at(array_index(get_id(elem), index_H)) = value;
            if (value > 0)
              values1_.at(array_index(get_id(elem), index_H)).resize(static_cast<std::size_t>(even_unknowns_on_node(static_cast<long>(this->get_expansion_order(elem, index_H)))));
            else
              values1_.at(array_index(get_id(elem), index_H)) = std::vector<ValueT>();
          }
          else if (topo_dim == flux_dim_)
          {
            expansion_order2_.at(array_index(get_id(elem), index_H)) = value;
            if (value > 0)
              values2_.at(array_index(get_id(elem), index_H)).resize(static_cast<std::size_t>(odd_unknowns_on_node(static_cast<long>(this->get_expansion_order(elem, index_H)))));
            else
              values2_.at(array_index(get_id(elem), index_H)) = std::vector<ValueT>();
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        std::size_t get_unknown_num(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            if ( this->get_unknown_index(elem, index_H) >= 0 )
              return static_cast<std::size_t>(even_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
            else
              return 0;
          }
          else if (topo_dim == flux_dim_)
          {
            if ( this->get_unknown_index(elem, index_H) >= 0 )
              return static_cast<std::size_t>(odd_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
            else
              return 0;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        bool   get_expansion_adaption(viennagrid_element_id elem) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return expansion_order_adaption_.at(get_id(elem));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }
        void set_expansion_adaption(viennagrid_element_id elem, bool b)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            expansion_order_adaption_.at(get_id(elem)) = b;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }



        // total energy
        ValueT get_value_H(std::size_t index_H) const { return values_H_.at(index_H); }
        void   set_value_H(std::size_t index_H, ValueT value) { values_H_.at(index_H) = value; }

        std::size_t get_value_H_size() const { return values_H_.size(); }

        // band gap center
        ValueT get_bandedge_shift(viennagrid_element_id elem) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return bandedge_shift1_.at(get_id(elem));
          }
          else if (topo_dim == flux_dim_)
          {
            return bandedge_shift2_.at(get_id(elem));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        void set_bandedge_shift(viennagrid_element_id elem, ValueT value)
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            bandedge_shift1_.at(get_id(elem)) = value;
          }
          else if (topo_dim == flux_dim_)
          {
            bandedge_shift2_.at(get_id(elem)) = value;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        ValueT get_kinetic_energy(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return (ctype_ == ELECTRON_TYPE_ID) ? get_value_H(index_H) - bandedge_shift1_.at(get_id(elem))
                                                : bandedge_shift1_.at(get_id(elem)) - get_value_H(index_H);
          }
          else if (topo_dim == flux_dim_)
          {
            return (ctype_ == ELECTRON_TYPE_ID) ? get_value_H(index_H) - bandedge_shift2_.at(get_id(elem))
                                                : bandedge_shift2_.at(get_id(elem)) - get_value_H(index_H);
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

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

        viennagrid_mesh mesh() const { return mesh_; }
        viennagrid_dimension density_dim() const { return density_dim_; }
        viennagrid_dimension    flux_dim() const { return flux_dim_; }

      private:
        std::size_t array_index(std::size_t element_id, std::size_t index_H) const { return element_id * values_H_.size() + index_H; }

        viennagrid_mesh                 mesh_;
        viennagrid_dimension            density_dim_;
        viennagrid_dimension            flux_dim_;
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
    template<typename ValueT = double>
    class const_she_quantity
    {
      std::size_t get_id(viennagrid_element_id elem) const { return std::size_t(viennagrid_index_from_element_id(elem)); }

      public:
        typedef ValueT         value_type;

        const_she_quantity(unknown_she_quantity<ValueT> const & quan)
          : mesh_(quan.mesh()),
            density_dim_(quan.density_dim()),
            flux_dim_(quan.flux_dim()),
            name_(quan.get_name()),
            values1_(quan.values1()),
            values2_(quan.values2()),
            expansion_order1_(quan.expansion_orders1()),
            expansion_order2_(quan.expansion_orders2()),
            get_value_H_size_(quan.get_value_H_size())
        {}

        ValueT const * get_values(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
            return &(values1_.at(array_index(get_id(elem), index_H)).at(0));
          else if (topo_dim == flux_dim_)
            return &(values2_.at(array_index(get_id(elem), index_H)).at(0));
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_values()");
        }

        std::size_t get_expansion_order(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            return expansion_order1_.at(array_index(get_id(elem), index_H));
          }
          else if (topo_dim == flux_dim_)
          {
            return expansion_order2_.at(array_index(get_id(elem), index_H));
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        std::size_t get_unknown_num(viennagrid_element_id elem, std::size_t index_H) const
        {
          viennagrid_dimension topo_dim = viennagrid_topological_dimension_from_element_id(elem);

          if (topo_dim == density_dim_)
          {
            if ( this->get_unknown_index(elem, index_H) >= 0 )
              return static_cast<std::size_t>(even_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
            else
              return 0;
          }
          else if (topo_dim == flux_dim_)
          {
            if ( this->get_unknown_index(elem, index_H) >= 0 )
              return static_cast<std::size_t>(odd_unknowns_on_node(static_cast<long>(get_expansion_order(elem, index_H))));
            else
              return 0;
          }
          else
            throw std::runtime_error("Invalid topological dimension for unknown_she_quantity::get_boundary_type()");
        }

        //ValueT operator()(AssociatedT1 const & elem, std::size_t index_H) const { return values1_.at(get_id(elem)).at(index_H); }

        std::string name() const { return name_; }

      private:
        std::size_t array_index(std::size_t element_id, std::size_t index_H) const { return element_id * get_value_H_size_ + index_H; }

        viennagrid_mesh                   mesh_;
        viennagrid_dimension              density_dim_;
        viennagrid_dimension              flux_dim_;
        std::string                       name_;
        std::vector<std::vector<ValueT> > values1_;
        std::vector<std::vector<ValueT> > values2_;
        std::vector<std::vector<short> >  expansion_order1_;
        std::vector<std::vector<short> >  expansion_order2_;
        std::size_t                       get_value_H_size_;

    };

  } //namespace she
} //namespace viennashe

#endif
