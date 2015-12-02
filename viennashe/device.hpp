#ifndef VIENNASHE_DEVICE_HPP
#define VIENNASHE_DEVICE_HPP
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

#include <stdexcept>

#include <ostream>
#include <cstddef>
#include <cmath>
#include <deque>

#include "viennagrid/viennagrid.h"

#include "viennashe/forwards.h"
#include "viennashe/exception.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/util/generate_device.hpp"

#include "viennashe/trap_level.hpp"

#include "viennashe/accessors.hpp"
#include "viennashe/setters.hpp"

/** @file  viennashe/device.hpp
    @brief Contains the definition of a device class independent of the actual macroscopic model to be solved.
*/

namespace viennashe
{

  /** @brief Defines the physical properties of a device, e.g. doping */
  template<typename MeshT>
  class device
  {
  public:
    typedef MeshT                           mesh_type;
    typedef viennagrid_region               segment_type;
    typedef viennagrid_region_id            segment_id_type;

    typedef long                            material_id_type;
    typedef std::size_t                     id_type;
    typedef trap_level                      trap_level_type;
    typedef std::vector<trap_level_type>    trap_level_container_type;

  private:
    std::size_t get_id(viennagrid_element_id el) const { return static_cast<std::size_t>(viennagrid_index_from_element_id(el)); }

  public:
    device() {}

    void load_mesh(std::string filename)
    {
      viennagrid_mesh_io mesh_io;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_create(&mesh_io));
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_mesh_set(mesh_io, mesh_));

      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_io_read(mesh_io, filename.c_str()));

      init_datastructures();
    }

    template<typename DeviceLoaderType>
    void load_device(DeviceLoaderType & loader)
    {
      loader(mesh_);
      init_datastructures();
    }


    void generate_mesh(viennashe::util::device_generation_config const & generator_params)
    {
      viennashe::util::generate_device(mesh_, generator_params);
      init_datastructures();
    }

    template<typename MeshGeneratorType>
    void generate_mesh(MeshGeneratorType const & gen)
    {
      gen(mesh_);
      init_datastructures();
    }

    /** @brief Returns the underlying mesh */
    MeshT const & mesh() const { return mesh_; }
    MeshT       & mesh()       { return mesh_; }

    segment_type segment(segment_id_type seg_id)
    {
      segment_type seg;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_region_get(mesh_, seg_id, &seg));
      return seg;
    }

    void scale(double factor)
    {
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_scale(mesh_, factor, NULL));
      init_datastructures(); //for Voronoi quantities
    }

    //
    // Lattice temperature: Cell-centric
    //

    /** @brief Sets the homogeneous temperature of the device. */
    void set_lattice_temperature(double new_value)
    {
      set_lattice_temp_region_impl(new_value, NULL);
    }

    /** @brief Sets the lattice temperature at a cell. */
    void set_lattice_temperature(double new_value, viennagrid_element_id c)
    {
      set_lattice_temp_impl(new_value, c);
    }

    /** @brief Sets the lattice temperature on a segment. */
    void set_lattice_temperature(double new_value, segment_type seg)
    {
      set_lattice_temp_region_impl(new_value, seg);
    }

    /** @brief Returns the lattice temperature on a cell */
    double get_lattice_temperature(viennagrid_element_id cell) const
    {
      assert(cell_temperature_.at(get_id(cell)) > 0 && bool("get_lattice_temp_impl(): Accessing non-positive temperature from cell!"));
      return cell_temperature_.at(get_id(cell));
    }

    /*double get_lattice_temperature(viennagrid_element_id facet) const
    {

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh_, facet, cell_dim, &cells_begin, &cells_end));

      if (cells_begin + 1 == cells_end)
        return get_lattice_temperature(*cells_begin);

      return (get_lattice_temperature(cells_begin[0]) + get_lattice_temperature(cells_begin[1])) / 2.0;
    }*/

    //
    // Doping: cell- and vertex-centric; cell-centric preferred. Conversion utilities are available.
    //

    /** @brief Sets the donator doping (in m^-3) in the specified cell */
    void set_doping_n(double value, viennagrid_element_id cell)
    {
      set_doping_n_impl(value, cell);
    }

    /** @brief Sets the donator doping (in m^-3) in the specified segment */
    void set_doping_n(double value, segment_type seg)
    {
      set_doping_np_impl(value, seg, DONATOR_DOPING_TYPE_ID);
    }

    /** @brief Sets the donator doping (in m^-3) in the whole device */
    void set_doping_n(double value)
    {
      set_doping_np_impl(value, NULL, DONATOR_DOPING_TYPE_ID);
    }

    /** @brief Returns the donator doping (in m^-3) in the specified cell */
    double get_doping_n(viennagrid_element_id cell_or_facet) const
    {
      return get_doping_element_np_impl(cell_or_facet, cell_doping_n_);
    }

    std::vector<double> const & doping_n() const { return cell_doping_n_; }


    /** @brief Sets the acceptor doping (in m^-3) in the specified cell */
    void set_doping_p(double value, viennagrid_element_id cell)
    {
      set_doping_p_impl(value, cell);
    }

    /** @brief Sets the acceptor doping (in m^-3) in the specified segment */
    void set_doping_p(double value, segment_type seg)
    {
      set_doping_np_impl(value, seg, ACCEPTOR_DOPING_TYPE_ID);
    }

    /** @brief Sets the acceptor doping (in m^-3) in the whole device */
    void set_doping_p(double value)
    {
      set_doping_np_impl(value, NULL, ACCEPTOR_DOPING_TYPE_ID);
    }

    /** @brief Returns the donator doping (in m^-3) in the specified cell */
    double get_doping_p(viennagrid_element_id cell_or_facet) const
    {
      return get_doping_element_np_impl(cell_or_facet, cell_doping_p_);
    }

    std::vector<double> const & doping_p() const { return cell_doping_p_; }

    double get_doping(viennagrid_element_id cell, carrier_type_id ctype) const
    {
      return (ctype == ELECTRON_TYPE_ID) ? get_doping_n(cell) : get_doping_p(cell);
    }

    //
    // Material: cell-centric.
    //
    /** @brief Sets the material ID on a cell */
    void set_material(long material_id, viennagrid_element_id elem)
    {
      cell_material_.at(get_id(elem)) = material_id;
    }

    /** @brief Sets the material type using the structs defined in viennashe::materials on a cell */
    template <typename MaterialType>
    void set_material(MaterialType, viennagrid_element_id elem)
    {
      set_material(long(MaterialType::id), elem);
    }

    //segment
    /** @brief Sets the material ID on a segment */
    void set_material(long material_id, segment_type seg)
    {
      set_material_on_complex(material_id, seg);
    }

    /** @brief Sets the material type using the structs defined in viennashe::materials on a segment */
    template <typename MaterialType>
    void set_material(MaterialType, segment_type seg)
    {
      set_material_on_complex(long(MaterialType::id), seg);
    }

    // full mesh
    /** @brief Sets a uniform material ID on the whole device */
    void set_material(long material_id)
    {
      set_material_on_complex(material_id, NULL);
    }

    /** @brief Sets a uniform material type using the structs defined in viennashe::materials on the whole device */
    template <typename MaterialType>
    void set_material(MaterialType)
    {
      set_material_on_complex(long(MaterialType::id), NULL);
    }

    /** @brief Returns the material id of the provided cell */
    long get_material(viennagrid_element_id cell) const
    {
      return cell_material_.at(get_id(cell));
    }

    std::vector<material_id_type> const & material() const { return cell_material_; }

    //
    // Contact potential: vertex-centric.
    //

  private:

    void init_datastructures()
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_int num_cells;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_element_count(mesh_, cell_dim, &num_cells));

      cell_doping_n_.resize(std::size_t(num_cells));
      cell_doping_p_.resize(std::size_t(num_cells));

      cell_temperature_.resize(std::size_t(num_cells), 300.0);

      cell_contact_potential_mask_.resize(std::size_t(num_cells));
      cell_contact_potential_.resize(std::size_t(num_cells), -1000.0);

      cell_material_.resize(std::size_t(num_cells));

      cell_traps_.resize(std::size_t(num_cells));

      cell_fixed_charges_.resize(std::size_t(num_cells));
    }

  public:
    /** @brief Sets the contact potential at a cell */
    void set_contact_potential(double pot, viennagrid_element_id cell)
    {
      cell_contact_potential_mask_.at(get_id(cell)) = true;
      cell_contact_potential_.at(get_id(cell)) = pot;
    }

    /** @brief Sets a contact potential for a whole segment */
    void set_contact_potential(double pot, segment_type seg)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh_, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        viennagrid_bool cell_in_segment;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(seg, *cit, &cell_in_segment));
        if (cell_in_segment)
          set_contact_potential(pot, *cit);
      }
    }

    /** @brief Returns the contact potential at a given cell (this is the externally applied voltage not considering any built-in potential) */
    double get_contact_potential(viennagrid_element_id cell) const
    {
      assert (cell_contact_potential_mask_.at(get_id(cell)) == true && bool("Accessing potential from vertex for which no contact potential was set!"));
      return cell_contact_potential_.at(get_id(cell));
    }

    /** @brief Returns true if a contact potential has been set for the respective vertex */
    bool has_contact_potential(viennagrid_element_id cell) const
    {
      return cell_contact_potential_mask_.at(get_id(cell));
    }

    //bool has_contact_potential(facet_type const &) const { return false; }

    //
    // Traps
    //

    /** @brief Adds a trap (density, energy) to a cell of the device */
    void add_trap_level(trap_level_type trap, viennagrid_element_id cell)
    {
      cell_traps_.at(get_id(cell)).push_back(trap);
    }

    /** @brief Adds a trap (density, energy) to a segment of the device */
    /*void add_trap_level(trap_level_type trap, segment_type const & seg)
    {
      add_trap_level_on_complex(trap, seg);
    }*/

    /** @brief Adds a trap (density, energy) to the whole device */
    void add_trap_level(trap_level_type trap)
    {
      add_trap_level_on_complex(trap, mesh_);
    }

    /** @brief Returns all the trap levels defined for the provided cell */
    trap_level_container_type const & get_trap_levels(viennagrid_element_id cell) const
    {
      return cell_traps_.at(get_id(cell));
    }

    /** @brief Removes all traps from the device */
    void clear_traps()
    {
      for (std::size_t i=0; i<cell_traps_.size(); ++i)
        cell_traps_[i].clear();
    }

    //
    // Fixed charges
    //

    /** @brief Sets a fixed charge at a cell.
     * @param charge The charge on the given cell. Use SI-unit: Coulomb
     * @param c The cell on which to put a fixed charge
     */
    void set_fixed_charge(viennagrid_element_id cell, double charge)
    {
      cell_fixed_charges_.at(get_id(cell)) = charge;
    }

    /**
     * @brief Gives the fixed charge set at a certain cell
     * @param c The cell
     * @return The charge in Coulomb
     */
    double get_fixed_charge(viennagrid_element_id cell) const
    {
      return cell_fixed_charges_.at(get_id(cell));
    }

protected:

    //
    // Lattice temperature
    //
    void set_lattice_temp_impl(double value, viennagrid_element_id cell)
    {
      cell_temperature_.at(get_id(cell)) = value;
    }

    /** @brief Sets the lattice temperature on either the full mesh or a particular segment.
     *
     *  @param value    Value of the lattice temperature in Kelvin
     *  @param mesh     The mesh
     *  @param seg      The segment (region). If NULL, then the value will be applied to the whole mesh.
     */
    void set_lattice_temp_region_impl(double value, segment_type seg)
    {
      if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_lattice_temp*: Lattice temperatures have to be greater 0 K!", value); }

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh_, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        viennagrid_bool cell_in_segment = true;
        if (seg)
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(seg, *cit, &cell_in_segment));
        if (cell_in_segment)
          set_lattice_temp_impl(value, *cit);
      }
    }


    //
    // Doping
    //

    void set_doping_n_impl(double value, viennagrid_element_id c)
    {
      if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping_p_impl(): Concentrations have to be greater 0 !", value); }
      cell_doping_n_.at(get_id(c)) = value;
    }

    void set_doping_p_impl(double value, viennagrid_element_id c)
    {
      if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping_p_impl(): Concentrations have to be greater 0 !", value); }
      cell_doping_p_.at(get_id(c)) = value;
    }

    void set_doping_np_impl(double value, segment_type seg, doping_type_id doping_type)
    {
      if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping*: Concentrations have to be greater 0 !", value); }

      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh_, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        viennagrid_bool cell_in_segment = true;
        if (seg)
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(seg, *cit, &cell_in_segment));
        if (cell_in_segment)
        {
          if (doping_type == DONATOR_DOPING_TYPE_ID)
            cell_doping_n_.at(get_id(*cit)) = value;
          else
            cell_doping_p_.at(get_id(*cit)) = value;
        }
      }
    }

    double get_doping_element_np_impl(viennagrid_element_id element, std::vector<double> const & doping_container) const
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      if (viennagrid_topological_dimension_from_element_id(element) == cell_dim) //cell
      {
        return doping_container.at(get_id(element));
      }
      else if (viennagrid_topological_dimension_from_element_id(element) == cell_dim - 1) // facet
      {
        viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
        VIENNASHE_VIENNAGRID_CHECK(viennagrid_element_coboundary_elements(mesh_, element, cell_dim, &cells_on_facet_begin, &cells_on_facet_end));

        if (cells_on_facet_begin + 1 == cells_on_facet_end)
          return doping_container.at(get_id(cells_on_facet_begin[0]));

        viennagrid_element_id c1 = cells_on_facet_begin[0];
        viennagrid_element_id c2 = cells_on_facet_begin[1];
        if (viennashe::materials::is_semiconductor(get_material(c1)))
        {
          if (viennashe::materials::is_semiconductor(get_material(c2)))
            return std::sqrt(doping_container.at(get_id(c1))) * std::sqrt(doping_container.at(get_id(c2)));
          else
            return doping_container.at(get_id(c1));
        }
        else if (viennashe::materials::is_semiconductor(get_material(c2)))
          return doping_container.at(get_id(c2));
        else
          return 0;
      }
      else
        throw std::runtime_error("get_doping_element_np_impl(): element is not cell or facet");
    }

    //
    // Material
    //

    void set_material_on_complex(long material_id, segment_type seg)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh_, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        viennagrid_bool cell_in_segment = true;
        if (seg)
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(seg, *cit, &cell_in_segment));
        if (cell_in_segment)
          set_material(material_id, *cit);
      }
    }



    //
    // Traps
    //

    void add_trap_level_on_complex(trap_level_type trap, segment_type seg)
    {
      viennagrid_dimension cell_dim;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_cell_dimension_get(mesh_, &cell_dim));

      viennagrid_element_id *cells_begin, *cells_end;
      VIENNASHE_VIENNAGRID_CHECK(viennagrid_mesh_elements_get(mesh_, cell_dim, &cells_begin, &cells_end));
      for (viennagrid_element_id *cit  = cells_begin;
                                  cit != cells_end;
                                ++cit)
      {
        viennagrid_bool cell_in_segment = true;
        if (seg)
          VIENNASHE_VIENNAGRID_CHECK(viennagrid_region_contains_element(seg, *cit, &cell_in_segment));
        if (cell_in_segment)
          add_trap_level(trap, *cit);
      }
    }

    MeshT             mesh_;

    // device data:
    std::vector<double>   cell_doping_n_;
    std::vector<double>   cell_doping_p_;

    std::vector<double>   cell_temperature_; //TODO: Think about handling temperatures in device

    std::vector<bool>     cell_contact_potential_mask_;
    std::vector<double>   cell_contact_potential_;

    std::vector<material_id_type> cell_material_;

    std::vector<trap_level_container_type> cell_traps_;

    std::vector<double>        cell_fixed_charges_;
  };


} //namespace viennashe

#endif
