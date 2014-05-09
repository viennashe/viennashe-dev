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

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/norm.hpp"
#include "viennagrid/algorithm/geometric_transform.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_reader.hpp"

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

  namespace detail
  {
    /** @brief Defines the physical properties of a device, e.g. doping */
    template <typename MeshT>
    class device_base
    {
      protected:
        typedef typename viennagrid::result_of::segmentation<MeshT>::type  segmentation_type;

        typedef typename viennagrid::result_of::vertex<MeshT>::type        vertex_type;
        typedef typename viennagrid::result_of::point<MeshT>::type         point_type;

        typedef typename viennagrid::result_of::const_cell_handle<MeshT>::type                       const_cell_handle_type;

      public:
        typedef MeshT                           mesh_type;
        typedef typename viennagrid::result_of::facet<MeshT>::type         facet_type;
        typedef typename viennagrid::result_of::cell<MeshT>::type          cell_type;

        typedef typename viennagrid::result_of::segment_handle<segmentation_type>::type  segment_type;
        typedef typename viennagrid::result_of::segmentation_segment_id_type<segmentation_type>::type   segment_id_type;
        typedef long                            material_id_type;
        typedef std::size_t                     id_type;
        typedef trap_level                      trap_level_type;
        typedef std::vector<trap_level_type>    trap_level_container_type;

        typedef typename viennagrid::result_of::voronoi_cell_contribution<const_cell_handle_type>::type   voronoi_contribution_container_type;

      private:
        std::size_t get_id(cell_type const & cell) const { return static_cast<std::size_t>(cell.id().get()); }

      public:
        device_base()
          : seg_(mesh_) {}

        void load_mesh(std::string filename)
        {
          if (filename.size() > 5 && filename.substr(filename.size()-5) == ".mesh")
          {
            viennagrid::io::netgen_reader my_reader;
            my_reader(mesh_, seg_, filename);
          }
          else if (filename.size() > 4 && (filename.substr(filename.size()-4) == ".vtu" || filename.substr(filename.size()-4) == ".pvd"))
          {
            viennagrid::io::vtk_reader<MeshT> mesh_reader;
            mesh_reader(mesh_, seg_, filename);
          }
          else
          {
            throw std::runtime_error("Unknown file extension!");
          }

          init_datastructures();
        }

        template <typename DeviceLoaderType>
        void load_device(DeviceLoaderType & loader)
        {
          loader(mesh_, seg_);
          init_datastructures();
        }


        void generate_mesh(viennashe::util::device_generation_config const & generator_params)
        {
          viennashe::util::generate_device(mesh_, seg_, generator_params);
          init_datastructures();
        }

        template < typename MeshGeneratorType >
        void generate_mesh(MeshGeneratorType const & gen)
        {
          gen(mesh_, seg_);
          init_datastructures();
        }

        /** @brief Returns the underlying mesh */
        MeshT const & mesh() const { return mesh_; }
        MeshT       & mesh()       { return mesh_; }

        segmentation_type const & segmentation() const { return seg_; }
        segmentation_type       & segmentation()       { return seg_; }

        segment_type const & segment(segment_id_type id) const { return seg_.at(id); }

        void scale(double factor)
        {
          viennagrid::scale(mesh_, factor);
          init_datastructures(); //for Voronoi quantities
        }

        //
        // Lattice temperature: Cell-centric
        //

        /** @brief Sets the homogeneous temperature of the device. */
        void set_lattice_temperature(double new_value)
        {
          set_lattice_temp_impl(new_value, mesh_);
        }

        /** @brief Sets the lattice temperature at a cell. */
        void set_lattice_temperature(double new_value, cell_type const & c)
        {
          set_lattice_temp_impl(new_value, c);
        }

        /** @brief Sets the lattice temperature on a segment. */
        void set_lattice_temperature(double new_value, segment_type const & s)
        {
          set_lattice_temp_impl(new_value, s);
        }

        /** @brief Returns the lattice temperature on a cell */
        double get_lattice_temperature(cell_type const & c) const
        {
          assert(cell_temperature_.at(get_id(c)) > 0 && bool("get_lattice_temp_impl(): Accessing non-positive temperature from cell!"));
          return cell_temperature_.at(get_id(c));
        }

        double get_lattice_temperature(facet_type const & facet) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<MeshT, facet_type, cell_type>::type     CellOnFacetContainer;
          typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

          CellOnFacetContainer cells_on_facet(mesh_, viennagrid::handle(mesh_, facet));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          cell_type const & c1 = *cofit;
          ++cofit;

          if (cofit == cells_on_facet.end())
            return get_lattice_temperature(c1);

          cell_type const & c2 = *cofit;
          return (get_lattice_temperature(c1) + get_lattice_temperature(c2)) / 2.0;
        }

        //
        // Doping: cell- and vertex-centric; cell-centric preferred. Conversion utilities are available.
        //

        /** @brief Sets the donator doping (in m^-3) in the specified cell */
        void set_doping_n(double value, cell_type const & c)
        {
          set_doping_n_impl(value, c);
        }

        /** @brief Sets the donator doping (in m^-3) in the specified segment */
        void set_doping_n(double value, segment_type const & d)
        {
          set_doping_n_impl(value, d);
        }

        /** @brief Sets the donator doping (in m^-3) in the specified segment */
        void set_doping_n(double value, segment_id_type const & seg_id)
        {
          set_doping_n_impl(value, segment(seg_id));
        }

        /** @brief Sets the donator doping (in m^-3) in the whole device */
        void set_doping_n(double value)
        {
          set_doping_n_impl(value, mesh_);
        }

        /** @brief Returns the donator doping (in m^-3) in the specified cell */
        double get_doping_n(cell_type const & c) const
        {
          return cell_doping_n_.at(get_id(c));
        }

        double get_doping_n(facet_type const & facet) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<MeshT, facet_type, cell_type>::type     CellOnFacetContainer;
          typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

          CellOnFacetContainer cells_on_facet(mesh_, viennagrid::handle(mesh_, facet));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          cell_type const & c1 = *cofit;
          ++cofit;

          if (cofit == cells_on_facet.end())
            return get_doping_n(c1);

          cell_type const & c2 = *cofit;
          if (viennashe::materials::is_semiconductor(get_material(c1)))
          {
            if (viennashe::materials::is_semiconductor(get_material(c2)))
              return std::sqrt(get_doping_n(c1)) * std::sqrt(get_doping_n(c2));
            else
              return get_doping_n(c1);
          }
          else if (viennashe::materials::is_semiconductor(get_material(c2)))
            return get_doping_n(c2);

          return 0;
        }

        std::vector<double> const & doping_n() const { return cell_doping_n_; }


        /** @brief Sets the acceptor doping (in m^-3) in the specified cell */
        void set_doping_p(double value, cell_type const & c)
        {
          set_doping_p_impl(value, c);
        }

        /** @brief Sets the acceptor doping (in m^-3) in the specified segment */
        void set_doping_p(double value, segment_type const & d)
        {
          set_doping_p_impl(value, d);
        }

        /** @brief Sets the donator doping (in m^-3) in the specified segment */
        void set_doping_p(double value, segment_id_type const & seg_id)
        {
          set_doping_p_impl(value, segment(seg_id));
        }

        /** @brief Sets the acceptor doping (in m^-3) in the whole device */
        void set_doping_p(double value)
        {
          set_doping_p_impl(value, mesh_);
        }

        /** @brief Returns the donator doping (in m^-3) in the specified cell */
        double get_doping_p(cell_type const & c) const
        {
          return cell_doping_p_.at(get_id(c));
        }

        double get_doping_p(facet_type const & facet) const
        {
          typedef typename viennagrid::result_of::const_coboundary_range<MeshT, facet_type, cell_type>::type     CellOnFacetContainer;
          typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

          CellOnFacetContainer cells_on_facet(mesh_, viennagrid::handle(mesh_, facet));

          CellOnFacetIterator cofit = cells_on_facet.begin();
          cell_type const & c1 = *cofit;
          ++cofit;

          if (cofit == cells_on_facet.end())
            return get_doping_p(c1);

          cell_type const & c2 = *cofit;
          if (viennashe::materials::is_semiconductor(get_material(c1)))
          {
            if (viennashe::materials::is_semiconductor(get_material(c2)))
              return std::sqrt(get_doping_p(c1)) * std::sqrt(get_doping_p(c2));
            else
              return get_doping_p(c1);
          }
          else if (viennashe::materials::is_semiconductor(get_material(c2)))
            return get_doping_p(c2);

          return 0;
        }

        std::vector<double> const & doping_p() const { return cell_doping_p_; }

        double get_doping(cell_type const & c, carrier_type_id ctype) const
        {
          return (ctype == ELECTRON_TYPE_ID) ? get_doping_n(c) : get_doping_p(c);
        }

        //
        // Material: cell-centric.
        //
        /** @brief Sets the material ID on a cell */
        void set_material(long material_id, cell_type const & elem)
        {
          cell_material_.at(get_id(elem)) = material_id;
        }

        /** @brief Sets the material type using the structs defined in viennashe::materials on a cell */
        template <typename MaterialType>
        void set_material(MaterialType, cell_type const & elem)
        {
          set_material(long(MaterialType::id), elem);
        }

        //segment
        /** @brief Sets the material ID on a segment */
        void set_material(long material_id, segment_type const & seg)
        {
          set_material_on_complex(material_id, seg);
        }

        /** @brief Sets the material ID on a segment */
        void set_material(long material_id, segment_id_type id)
        {
          set_material_on_complex(material_id, segment(id));
        }

        /** @brief Sets the material type using the structs defined in viennashe::materials on a segment */
        template <typename MaterialType>
        void set_material(MaterialType, segment_type const & seg)
        {
          set_material_on_complex(long(MaterialType::id), seg);
        }

        // full mesh
        /** @brief Sets a uniform material ID on the whole device */
        void set_material(long material_id)
        {
          set_material_on_complex(material_id, mesh_);
        }

        /** @brief Sets a uniform material type using the structs defined in viennashe::materials on the whole device */
        template <typename MaterialType>
        void set_material(MaterialType)
        {
          set_material_on_complex(long(MaterialType::id), mesh_);
        }

        /** @brief Returns the material id of the provided cell */
        long get_material(cell_type const & elem) const
        {
          return cell_material_.at(get_id(elem));
        }

        std::vector<material_id_type> const & material() const { return cell_material_; }

        //
        // Contact potential: vertex-centric.
        //

      private:

        void init_datastructures()
        {
          cell_doping_n_.resize(viennagrid::cells(mesh_).size());
          cell_doping_p_.resize(viennagrid::cells(mesh_).size());

          cell_temperature_.resize(viennagrid::cells(mesh_).size(), 300.0);

          cell_contact_potential_mask_.resize(viennagrid::cells(mesh_).size());
          cell_contact_potential_.resize(viennagrid::cells(mesh_).size(), -1000.0);

          cell_material_.resize(viennagrid::cells(mesh_).size());

          cell_traps_.resize(viennagrid::cells(mesh_).size());

          cell_fixed_charges_.resize(viennagrid::cells(mesh_).size());
        }

      public:
        /** @brief Sets the contact potential at a cell */
        void set_contact_potential(double pot, cell_type const & c)
        {
          cell_contact_potential_mask_.at(get_id(c)) = true;
          cell_contact_potential_.at(get_id(c)) = pot;
        }

        /** @brief Sets a contact potential for a whole segment */
        void set_contact_potential(double pot, segment_type const & seg)
        {
          typedef typename viennagrid::result_of::const_cell_range<segment_type>::type        CellOnSegmentContainer;
          typedef typename viennagrid::result_of::iterator<CellOnSegmentContainer>::type      CellOnSegmentIterator;

          CellOnSegmentContainer cells_on_segment(seg);
          for (CellOnSegmentIterator cit  = cells_on_segment.begin();
                                     cit != cells_on_segment.end();
                                   ++cit)
          {
            set_contact_potential(pot, *cit);
          }
        }

        /** @brief Returns the contact potential at a given cell (this is the externally applied voltage not considering any built-in potential) */
        double get_contact_potential(cell_type const & c) const
        {
          assert (cell_contact_potential_mask_.at(get_id(c)) == true && bool("Accessing potential from vertex for which no contact potential was set!"));
          return cell_contact_potential_.at(get_id(c));
        }

        /** @brief Returns true if a contact potential has been set for the respective vertex */
        bool has_contact_potential(cell_type const & c) const
        {
          return cell_contact_potential_mask_.at(get_id(c));
        }

        bool has_contact_potential(facet_type const &) const { return false; }

        //
        // Traps
        //

        /** @brief Adds a trap (density, energy) to a cell of the device */
        void add_trap_level(trap_level_type trap, cell_type const & cell)
        {
          cell_traps_.at(get_id(cell)).push_back(trap);
        }

        /** @brief Adds a trap (density, energy) to a segment of the device */
        void add_trap_level(trap_level_type trap, segment_type const & seg)
        {
          add_trap_level_on_complex(trap, seg);
        }

        /** @brief Adds a trap (density, energy) to the whole device */
        void add_trap_level(trap_level_type trap)
        {
          add_trap_level_on_complex(trap, mesh_);
        }

        /** @brief Returns all the trap levels defined for the provided cell */
        trap_level_container_type const & get_trap_levels(cell_type const & cell) const
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
        void set_fixed_charge(cell_type const & c, double charge)
        {
          cell_fixed_charges_.at(get_id(c)) = charge;
        }

        /**
         * @brief Gives the fixed charge set at a certain cell
         * @param c The cell
         * @return The charge in Coulomb
         */
        double get_fixed_charge(cell_type const & c) const
        {
          return cell_fixed_charges_.at(get_id(c));
        }

    protected:

        //
        // Lattice temperature
        //
        void set_lattice_temp_impl(double value, cell_type const & c)
        {
          cell_temperature_.at(get_id(c)) = value;
        }

        template <typename MeshOrSegmentT>
        void set_lattice_temp_impl(double value, MeshOrSegmentT const & meshseg)    //ComplexType is a segment or a mesh
        {
          typedef typename viennagrid::result_of::const_cell_range<MeshOrSegmentT>::type     CellContainer;
          typedef typename viennagrid::result_of::iterator<CellContainer>::type              CellIterator;
          if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_lattice_temp*: Lattice temperatures have to be greater 0 K!", value); }

          CellContainer cells(meshseg);
          for (CellIterator cit  = cells.begin();
                            cit != cells.end();
                          ++cit)
          {
            set_lattice_temp_impl(value, *cit);
          }
        }


        //
        // Doping
        //

        void set_doping_n_impl(double value, cell_type const & c)
        {
          if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping_p_impl(): Concentrations have to be greater 0 !", value); }
          cell_doping_n_.at(get_id(c)) = value;
        }

        void set_doping_p_impl(double value, cell_type const & c)
        {
          if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping_p_impl(): Concentrations have to be greater 0 !", value); }
          cell_doping_p_.at(get_id(c)) = value;
        }

        template <typename MeshOrSegmentT>
        void set_doping_n_impl(double value, MeshOrSegmentT const & meshseg)
        {
          typedef typename viennagrid::result_of::const_cell_range<MeshOrSegmentT>::type   CellContainer;
          typedef typename viennagrid::result_of::iterator<CellContainer>::type            CellIterator;

          if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping*: Concentrations have to be greater 0 !", value); }
          CellContainer cells(meshseg);
          for (CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit)
          {
            set_doping_n_impl(value, *cit);
          }
        }

        template <typename MeshOrSegmentT>
        void set_doping_p_impl(double value, MeshOrSegmentT const & meshseg)
        {
          typedef typename viennagrid::result_of::const_cell_range<MeshOrSegmentT>::type    CellContainer;
          typedef typename viennagrid::result_of::iterator<CellContainer>::type             CellIterator;

          if (value <= 0.0) { throw viennashe::invalid_value_exception("device.set_doping*: Concentrations have to be greater 0 !", value); }
          CellContainer cells(meshseg);
          for (CellIterator cit = cells.begin();
              cit != cells.end();
              ++cit)
          {
            set_doping_p_impl(value, *cit);
          }
        }

        //
        // Material
        //

        template <typename MeshOrSegmentT>
        void set_material_on_complex(long material_id, MeshOrSegmentT const & meshseg)
        {
          typedef typename viennagrid::result_of::const_cell_range<MeshOrSegmentT>::type    CellContainer;
          typedef typename viennagrid::result_of::iterator<CellContainer>::type             CellIterator;
          //if (material_id < 0.0) { throw viennashe::invalid_value_exception("device.set_material*:  !", material_id); }

          CellContainer cells(meshseg);
          for (CellIterator cit = cells.begin();
                            cit != cells.end();
                          ++cit)
            set_material(material_id, *cit);
        }



        //
        // Traps
        //

        template <typename MeshOrSegmentT>
        void add_trap_level_on_complex(trap_level_type trap, MeshOrSegmentT const & meshseg)
        {
          typedef typename viennagrid::result_of::const_cell_range<MeshOrSegmentT>::type   CellContainer;
          typedef typename viennagrid::result_of::iterator<CellContainer>::type               CellIterator;

          CellContainer cells(meshseg);
          for (CellIterator cit = cells.begin();
                            cit != cells.end();
                          ++cit)
          {
            add_trap_level(trap, *cit);
          }
        }

        MeshT             mesh_;
        segmentation_type seg_;

        // device data:
        std::vector<double>   cell_doping_n_;
        std::vector<double>   cell_doping_p_;

        std::vector<double>         cell_temperature_; //TODO: Think about handling temperatures in device

        std::vector<bool>     cell_contact_potential_mask_;
        std::vector<double>   cell_contact_potential_;

        std::vector<material_id_type> cell_material_;

        std::vector<trap_level_container_type> cell_traps_;

        std::vector<double>        cell_fixed_charges_;
    };

  }


  /** @brief Defines the physical properties of a device, e.g. doping. This is the implementation for 2d and higher dimensions */
  template <typename MeshT,
            bool edges_and_cells_different /* see forwards.h for default argument */ >
  class device : public detail::device_base<MeshT> {};



} //namespace viennashe

#endif
