#ifndef VIENNASHE_IO_ADD_TO_WRITER_HPP
#define	VIENNASHE_IO_ADD_TO_WRITER_HPP
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

/** @file viennashe/io/add_to_writer.hpp
    @brief Convenience functions for viennagrid::io (mostly adders)
*/


// std
#include <string>
#include <vector>

// viennashe
#include "viennashe/forwards.h"

#include "viennashe/postproc/current_density.hpp"
#include "viennashe/postproc/electric_field.hpp"
#include "viennashe/she/postproc/all.hpp"


namespace viennashe
{
  namespace io
  {
    namespace detail
    {

      /**
       * @brief A wrapper for the accessor in viennashe, in order for them to be suitable for viennagrid.
       * @warning The ViennaSHE accessors are held as <b>copies</b>!
       */
      template < typename AccessorT, typename ValueT = typename AccessorT::value_type >
      class viennagrid_accessor_wrapper
      {
      public:
        typedef ValueT value_type;

        viennagrid_accessor_wrapper(AccessorT acc) : acc_(acc) { }
        viennagrid_accessor_wrapper(viennagrid_accessor_wrapper const & o) : acc_(o.acc_) { }

        template < typename ElementT >
        value_type const & at(ElementT const & el) const
        {
          static value_type cache_ = acc_(el);
          cache_ = acc_(el);
          return cache_;
        }

        template < typename ElementT >
        value_type const & operator()(ElementT const & el) const
        {
          static value_type cache_ = acc_(el);
          cache_ = acc_(el);
          return cache_;
        }

        template < typename ElementT >
        value_type const * find(ElementT const & el) const
        {
          static value_type cache_ = acc_(el);
          cache_ = acc_(el);
          return &cache_;
        }

        private:
          AccessorT acc_; // COPY !!!
      };

      /**
       * @brief Implementation of the quantity addition to a writer for single valued quantities
       * @param quantity An accessor to the quantity
       * @param writer The writer, e.g. vtk writer
       * @param name The name of the quantity
       * @param dummy Do not touch! This is a tag.
       */
      template <typename MacroscopicQuantityAccessor, typename WriterType>
      void add_macroscopic_quantity_to_writer_impl(MacroscopicQuantityAccessor const & quantity,
                                              WriterType & writer,
                                              std::string name,
                                              double * dummy = 0)
      {
        (void)dummy;
        // quantity IS GOING TO BE COPIED!
        viennagrid_accessor_wrapper<MacroscopicQuantityAccessor, double> wrap(quantity);

        writer.add_scalar_data_on_cells(wrap, name);
      }

      /**
       * @brief Implementation of the quantity addition to a writer for vector valued quantities
       * @param quantity An accessor to the quantity
       * @param writer The writer, e.g. vtk writer
       * @param name The name of the quantity
       * @param dummy Do not touch! This is a tag.
       */
      template <typename MacroscopicQuantityAccessor, typename WriterType>
      void add_macroscopic_quantity_to_writer_impl(MacroscopicQuantityAccessor const & quantity,
                                              WriterType & writer,
                                              std::string name,
                                              std::vector<double> * dummy = 0)
      {
        (void)dummy;
        // quantity IS GOING TO BE COPIED!
        viennagrid_accessor_wrapper<MacroscopicQuantityAccessor, std::vector<double> > wrap(quantity);
        writer.add_vector_data_on_cells(wrap, name);
      }
    } // namespace detail

    /**
     * @brief Adds a quantity to a writer
     * @param quantity The accessor to the quantity
     * @param writer The writer
     * @param name Name of the quantity in the file
     */
    template <typename MacroscopicQuantityAccessor, typename WriterType>
    void add_macroscopic_quantity_to_writer(MacroscopicQuantityAccessor const & quantity,
                                            WriterType & writer, std::string name)
    {
      typename MacroscopicQuantityAccessor::value_type * dummy = 0;
       // quantity IS GOING TO BE COPIED!
      viennashe::io::detail::add_macroscopic_quantity_to_writer_impl(quantity, writer, name, dummy);
    }

    /**
     * @brief Special free function to add the potential to a writer
     * @param pot The accessor to the electrostatic potential
     * @param writer The writer
     * @param name Name of the elec. potential in the file
     */
    template <typename PotentialAccessor,
              typename WriterType>
    void add_potential_to_writer(PotentialAccessor const & pot,
                                 WriterType & writer, std::string name)
    {
       // quantity IS GOING TO BE COPIED!
      viennashe::io::add_macroscopic_quantity_to_writer(pot, writer, name);
    }

    /*
    template <typename DeviceType,
              typename SHEQuantity,
              typename WriterType>
    void add_carrier_velocity_to_writer(DeviceType  const & device,
                                        SHEQuantity const & quan,
                                        viennashe::config::dispersion_relation_type const & dispersion,
                                        WriterType & writer,
                                        std::string name)
    {
      // get accessor
      viennashe::she::carrier_velocity_wrapper<DeviceType, SHEQuantity> velocity_wrapper(device, quan, dispersion);
      // write
      // quantity IS GOING TO BE COPIED!
      viennashe::io::add_macroscopic_quantity_to_writer(velocity_wrapper, writer, name);
    } */

    /**
     * @brief Special free function to add the current density - calculated from a solution of the BTE - to a writer
     * @param device The device
     * @param quan   The SHE quantity, i.e. the EDF
     * @param conf   The simulator configuration
     * @param writer The writer
     * @param name   The name of the current density in the file
     */
    template <typename SHEQuantity,
              typename WriterType>
    void add_current_density_to_writer(viennashe::device const & device,
                                       viennashe::config const & conf,
                                       SHEQuantity const & quan,
                                       WriterType & writer,
                                       std::string name)
    {
      viennashe::she::current_density_wrapper<SHEQuantity> current_wrapper(device, conf, quan);
      // quantity IS GOING TO BE COPIED!
      viennashe::io::add_macroscopic_quantity_to_writer(current_wrapper, writer, name);
    }

    /**
     * @brief Special free function to add the current density - calculated from DD - to a writer
     * @param device The device
     * @param potential Accessor to the elec. potential
     * @param carrier Accessor to the carrier density
     * @param ctype Carrier type (electrons or holes)
     * @param mobility_model The mobility model
     * @param writer The writer
     * @param name The name of the current density in the file
    */
    template <typename PotentialQuantityType,
              typename CarrierQuantityType,
              typename MobilityModel,
              typename WriterType>
    void add_current_density_to_writer(viennashe::device const & device,
                                       PotentialQuantityType const & potential,
                                       CarrierQuantityType   const & carrier,
                                       viennashe::carrier_type_id ctype,
                                       MobilityModel const & mobility_model,
                                       WriterType & writer,
                                       std::string name)
    {
      typedef typename viennashe::current_density_wrapper<PotentialQuantityType,
                                                          CarrierQuantityType, MobilityModel>  current_density_type;
      current_density_type Jfield(device, ctype, potential, carrier, mobility_model);
      // quantity IS GOING TO BE COPIED!
      viennashe::io::add_macroscopic_quantity_to_writer(Jfield, writer, name);
    }

    /**
     * @brief Free function to add the average carrier energy to a writer
     * @param quan The SHE quantity, that is the EDF
     * @param conf The simulator configuration
     * @param writer The writer
     * @param name Name of the average carrier energy in the file
     */
    template <typename SHEQuantity,
              typename WriterType>
    void add_kinetic_carrier_energy_to_writer(viennashe::config const & conf,
                                              SHEQuantity const & quan,
                                              WriterType & writer,
                                              std::string name)
    {
      viennashe::she::carrier_energy_wrapper<SHEQuantity> wrapper(conf, quan);
      // quantity IS GOING TO BE COPIED!
      viennashe::io::add_macroscopic_quantity_to_writer(wrapper, writer, name);
    }


  } // namespace io
} // namespace viennashe


#endif	/* VIENNASHE_IO_ADD_TO_WRITER_HPP */

