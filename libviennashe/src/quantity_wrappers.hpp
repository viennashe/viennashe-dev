#ifndef LIBVIENNASHE_QUANTITY_WRAPPERS_HPP
#define	LIBVIENNASHE_QUANTITY_WRAPPERS_HPP
/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */


/** @file libviennashe/src/quantity_wrappers.hpp
    @brief Routines to wrap ViennaSHE quantities (doping, potential, distribution functions, ...)
*/

// ++++++++++++++++++++++++++++++++++++
//
// +++++++++++   C++ ONLY   +++++++++++
//
// ++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cstdlib>
#include <vector>

#include "viennashe/forwards.h"
#include "viennashe/exception.hpp"

#include "libviennashe/src/viennashe_all.hpp"


namespace libviennashe
{
  /** @brief Contains everything quantity related */
  namespace quantity
  {
    /** @brief The main quantity wrapper <b>interface</b> */
    class quantity_wrapper
    {
    public:
      quantity_wrapper(std::string name) : name_(name) {}

      /** @brief Default CTOR. Creates an empty wrapper. */
      quantity_wrapper() : name_("") {}

      quantity_wrapper(quantity_wrapper const & v) : name_(v.name_) {}
      virtual ~quantity_wrapper() {}

      void operator=(quantity_wrapper const & o) { this->name_ = o.name_; }

      /**
       * @brief Interface. Fills values with all quantity values,
       *        where the length of values[i] is to be found in len[i]
       * @param values The quantity values
       * @param len The length of each quantity value (1 = scalar)
       */
      virtual void fill(double ** values, viennashe_index_type * len) const = 0;

      /**
       * @brief Interface. Fills values with a single value (scalar or vector) at the given element index
       * @param idx The index of the element
       * @param values A std::vector<double> holding the values
       */
      virtual void fill_single (std::size_t idx, std::vector<double> & values) const  = 0;

      /** @brief A simple copy factory for storage */
      virtual quantity_wrapper * copy() const = 0;

      /** @brief Returns the unique name of the quantity */
      std::string name() const { return this->name_; }

    protected:
      /** @brief Name setter for implementing classes */
      void set_name(std::string const & n) { this->name_ = n; }

    private:
      std::string name_;
    };

    /** @brief Implements quantity_wrapper. Wraps scalar quantities, which are accessible via a
     *         device based accessor in ViennaSHE (cf. accessor.hpp) */
    template <typename DeviceT, typename AccessorT, typename ElementTagT >
    class accessor_based_quantity_wrapper : public quantity_wrapper
    {
    public:

      typedef accessor_based_quantity_wrapper<DeviceT, AccessorT, ElementTagT> self_type;

      /**
       * @brief CTOR.
       * @param acc The accessor. <b>Must be deep copy-able!</b>
       * @param dev The device. A const reference to the device will be stored by this class
       * @param name The unique name of the quantity
       */
      accessor_based_quantity_wrapper(AccessorT const & acc,
        DeviceT const & dev,
        std::string name)
        : quantity_wrapper(name), dev_(dev), acc_(acc)
      { }

      accessor_based_quantity_wrapper(accessor_based_quantity_wrapper const & o)
        : quantity_wrapper(o), dev_(o.dev_), acc_(o.acc_)
      { }

      void operator=(accessor_based_quantity_wrapper const & o)
      {
        this->set_name(o.name());
        dev_(o.dev_); acc_(o.acc_);
      }

      /** @brief Simple forward to the accessor.
       * @param id The id of the element on which the data is stored
       */
      double get(std::size_t id) const
      {
        return acc_(viennagrid::elements<ElementTagT>(dev_.mesh())[id]);
      }

      /**
       * @brief Implementation of fill.
       */
      virtual void fill(double ** values, viennashe_index_type * len) const
      {
        const std::size_t num = viennagrid::elements<ElementTagT>(dev_.mesh()).size();

        for (std::size_t i = 0; i < num; ++i)
        {
          values[i] = new double(); // The user has to clear things up
          values[i][0] = this->get(i);
          len[i] = 1;
        }
      }

      /**
       * @brief Implemenation of fill_single. Uses get()
       */
      virtual void fill_single(std::size_t idx, std::vector<double> & values) const
      {
        values.resize(1);
        values[0] = this->get(idx);
      }

      virtual quantity_wrapper * copy() const { return new self_type(*this); }

    private:
      DeviceT   const & dev_;
      AccessorT acc_;
    };

    /** @brief Implements quantity_wrapper. Wraps vector valued quantities, which are accessible via a
     *         device based accessor in ViennaSHE (cf. accessor.hpp)
     */
    template <typename DeviceT, typename AccessorT, typename ElementTagT >
    class accessor_based_array_quantity_wrapper : public quantity_wrapper
    {
    public:

      typedef accessor_based_array_quantity_wrapper<DeviceT, AccessorT, ElementTagT> self_type;

      /**
       * @brief CTOR.
       * @param acc The accessor. <b>Must be deep copy-able!</b>
       * @param dev The device. A reference to this object will be held
       * @param name The unique name of the vector valued quantity
       */
      accessor_based_array_quantity_wrapper(AccessorT const & acc,
        DeviceT const & dev,
        std::string name)
        : quantity_wrapper(name), dev_(dev), acc_(acc)
      { }

      accessor_based_array_quantity_wrapper(accessor_based_array_quantity_wrapper const & o)
        : quantity_wrapper(o), dev_(o.dev_), acc_(o.acc_)
      { }

      void operator=(accessor_based_array_quantity_wrapper const & o)
      {
        this->set_name(o.name());
        dev_(o.dev_); acc_(o.acc_);
      }

      /**
       * @brief Simple forward to the accessor
       * @param id The elment id of the element on which the vector is stored on
       * @return A vector as std::vector<double>
       */
      std::vector<double> get(std::size_t id) const
      {
        return acc_(viennagrid::elements<ElementTagT>(dev_.mesh())[id]);
      }

      /** @brief Implements fill_single. */
      virtual void fill_single(std::size_t idx, std::vector<double> & values) const { values = this->get(idx); }

      /** @brief Implements fill */
      virtual void fill(double ** values, viennashe_index_type * len) const
      {
        const std::size_t num = viennagrid::elements<ElementTagT>(dev_.mesh()).size();

        for (std::size_t i = 0; i < num; ++i)
        {
          const std::vector<double> valuevector = this->get(i);
          const std::size_t current_length = valuevector.size();

          values[i] = new double[current_length]; // The user has to clear things up

          for (std::size_t j = 0; j < current_length; ++j) values[i][j] = valuevector.at(j);

          len[i] = static_cast<viennashe_index_type>(current_length);
        }
      }

      /** @brief Implements copy */
      virtual quantity_wrapper * copy() const { return new self_type(*this); }

    private:
      DeviceT   const & dev_;
      AccessorT acc_;
    };

    /** @brief The main quantity register. Not deep copy-able! */
    class quantitiy_register
    {
    public:
      typedef std::map<std::string, quantity_wrapper * > quan_map;

      virtual ~quantitiy_register()
      {
        for (quan_map::iterator it = quans_.begin(); it != quans_.end(); ++it)
        {
         if (it->second != NULL)
           delete it->second;
         it->second = NULL;
        }
      }

      /** @brief Returns the number of registered quantities */
      std::size_t size() const { return this->quans_.size(); }

      /** @brief Returns all names of the registered quantities */
      std::vector<std::string> get_names() const
      {
        std::vector<std::string> v;
        for (quan_map::const_iterator it = quans_.begin(); it != quans_.end(); ++it)
        {
          v.push_back(it->first);
        }
        return v;
      }

      /** @brief Returns a single const reference to a quantitiy wrapper
       * @param name The name of the quantity to retrieve
       */
      quantity_wrapper const & get(std::string const & name) const { return *(quans_.at(name)); }

      /**
       * @brief Returns true if a quantity is registered
       * @param name The name of the quantity
       * @return true if the quantity exists (is registered)
       */
      bool has_quan(std::string const & name) const { return (this->quans_.find(name) != this->quans_.end()); }

      /** @brief Registers a quantity
       * @param quan The quantity by reference. Will be copied using quan.copy()!
       */
      void register_quan(quantity_wrapper const & quan)
      {
        std::string name = quan.name();
        if (this->has_quan(name))
        {
          viennashe::log::error() << "ERROR: ViennaSHE C/C++ Wrappers failed! Where? Quantitiy Management at register_quan for '"
                                  << name << "', which is already registered!" << std::endl;
          throw std::runtime_error("Quantitiy already exists");
        }
        else
        {
          this->quans_.insert(std::make_pair(name, quan.copy()));
        }
      }

    private:
      quan_map quans_;
    }; // quantity_register


  } // namespace quantity

  /** @brief C++ to C wrapper of the quantity registry */
  struct quan_register_internal
  {

    quan_register_internal() : int_sim(NULL) { }

    std::size_t count() const
    {
      return cell_based.size() ;
    }

    libviennashe::quantity::quantitiy_register cell_based; //! Register for vertex based quantities

    viennashe_simulator_impl * int_sim; //! Simulator. Not deleted by this wrapper
  };

} // namespace libviennashe

#endif /* LIBVIENNASHE_QUANTITY_WRAPPERS_HPP */

