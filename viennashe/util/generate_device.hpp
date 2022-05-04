#ifndef VIENNASHE_UTIL_GENERATE_DEVICE_HPP
#define VIENNASHE_UTIL_GENERATE_DEVICE_HPP

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

// std
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/mesh/element_creation.hpp"
#include "viennagrid/topology/line.hpp"
#include "viennagrid/topology/quadrilateral.hpp"

// viennashe
#include "viennashe/log/log.hpp"
#include "viennashe/util/log_keys.h"


/** @file viennashe/util/generate_device.hpp
    @brief Contains a very simple mesh generator (ortho-grids) for one and two spatial dimensions.
*/

namespace viennashe
{
  namespace util
  {


    /** @brief Configuration class for the simple mesh generator */
    class device_generation_config
    {
    public:
      class segment_description
      {
      public:
        segment_description()
          : start_x_(0.0),
            start_y_(0.0),
            len_x_(0.0),
            len_y_(0.0),
            points_x_(1),
            points_y_(1) {}

        segment_description(double start_x, double len_x, unsigned long points_x)
          : start_x_(start_x),
            start_y_(0.0),
            len_x_(len_x),
            len_y_(0.0),
            points_x_(points_x),
            points_y_(1) {}

        segment_description(double start_x, double start_y, double len_x, double len_y, unsigned long points_x, unsigned long points_y)
          : start_x_(start_x),
            start_y_(start_y),
            len_x_(len_x),
            len_y_(len_y),
            points_x_(points_x),
            points_y_(points_y) {}

        double get_start_x() const { return start_x_; }
        double get_start_y() const { return start_y_; }

        double get_length_x() const { return len_x_; }
        double get_length_y() const { return len_y_; }

        unsigned long get_points_x() const { return points_x_; }
        unsigned long get_points_y() const { return points_y_; }

      private:
        double start_x_;  /// Length in x-direction.
        double start_y_;  /// Length in y-direction. Ignored for one-dimensional devices.

        double len_x_;    /// Length in x-direction.
        double len_y_;    /// Length in y-direction. Zero for one-dimensional devices

        unsigned long points_x_;  /// Number of points in x-direction.
        unsigned long points_y_;  /// Number of points in y-direction. Zero for one-dimensional devices
      };

      typedef segment_description   segment_description_type;

      void clear() { segment_descs_.clear(); }

      void add_segment(double start_x, double start_y, double len_x, double len_y, unsigned long points_x, unsigned long points_y)
      {
        segment_descs_.push_back(segment_description(start_x, start_y, len_x, len_y, points_x, points_y));
      }

      void add_segment(double start_x, double len_x, unsigned long points_x,
                       double start_y, double len_y, unsigned long points_y)
      {
        segment_descs_.push_back(segment_description(start_x, start_y, len_x, len_y, points_x, points_y));
      }

      void add_segment(double start_x, double len_x, unsigned long points_x)
      {
        add_segment(start_x, 0.0, len_x, 0.0, points_x, 1);
      }

      std::size_t size() const { return segment_descs_.size(); }
      segment_description const & at(std::size_t i) const { return segment_descs_.at(i); }

    private:
      std::vector<segment_description>  segment_descs_;
    };

    namespace detail
    {
      //
      // 1d generation
      //

      /** @brief Implementation of a simple one-dimensional 'mesh' generation. Fills a ViennaGrid mesh.
       *
       * @param mesh           An empty ViennaGrid mesh
       * @param segmentation   A mesh segmentation for setting up the segments
       * @param conf           The mesh generation configuration object
       */
      template <typename MeshT, typename SegmentationT>
      void generate_device_impl(MeshT & mesh,
                                SegmentationT & segmentation,
                                device_generation_config const & conf,
                                viennagrid::simplex_tag<1>
                               )
      {
        typedef typename viennagrid::result_of::point<MeshT>::type     PointType;
        typedef typename viennagrid::result_of::vertex<MeshT>::type    VertexType;
        typedef typename viennagrid::result_of::cell<MeshT>::type      CellType;

        typedef typename viennagrid::result_of::cell_tag<MeshT>::type  CellTag;

        typedef typename viennagrid::result_of::handle<MeshT, viennagrid::vertex_tag>::type  VertexHandleType;

        typedef typename device_generation_config::segment_description   SegmentDescriptionType;


        //
        // Prepare vertices at segment boundaries (note that this is O(N^2) with respect to the number of segments N. Not expected to hurt in practice, though...):
        //
        std::vector<PointType> segment_boundary_points;

        std::size_t total_points = 1;
        for (std::size_t i=0; i<conf.size(); ++i)
        {
          SegmentDescriptionType const & seg_desc = conf.at(i);

          assert(seg_desc.get_points_x() > 1   && bool("Logic error: Not enough points in x-direction provided for segment"));
          assert(seg_desc.get_points_y() == 1  && bool("Logic error: Provided two-dimensional grid description for one-dimensional mesh"));
          assert(seg_desc.get_length_x() > 0.0 && bool("Logic error: x-coordinate is degenerate in device generation."));

          double distance_tolerance = 1e-10 * seg_desc.get_length_x();

          PointType p0(seg_desc.get_start_x());
          PointType p1(seg_desc.get_start_x() + seg_desc.get_length_x());

          bool insert_p0 = true;
          bool insert_p1 = true;

          for (std::size_t j=0; j<segment_boundary_points.size(); ++j)
          {
            if (viennagrid::norm_2(p0 - segment_boundary_points[j]) < distance_tolerance )
              insert_p0 = false;
            if (viennagrid::norm_2(p1 - segment_boundary_points[j]) < distance_tolerance )
              insert_p1 = false;
          }

          if (insert_p0)
            segment_boundary_points.push_back(p0);
          if (insert_p1)
            segment_boundary_points.push_back(p1);

          total_points += seg_desc.get_points_x() - 1;
        }

        //
        // Set up the vertices:
        //
        log::info<log_generate_device>() << "* generate_device(): Setting up vertices..." << std::endl;

        std::vector<int> segment_boundary_ids(segment_boundary_points.size(), -1);

        std::vector< std::vector<int> > segment_vertex_ids(conf.size()); // store the global vertex ID for each segment to save lookups later on.

        int vertex_counter = 0;

        // iterate over all segments and insert points if not already inserted in mesh. Store vertex ID for each segment.
        for (std::size_t i=0; i<conf.size(); ++i)
        {
          SegmentDescriptionType const & seg_desc = conf.at(i);

          segment_vertex_ids[i].resize(seg_desc.get_points_x());

          double distance_tolerance = 1e-10 * seg_desc.get_length_x();

          // Add get_points_x() points in the segment to the mesh if they haven't been added yet (check segment boundaries)
          for (std::size_t j = 0; j<seg_desc.get_points_x(); ++j)
          {
            PointType candidate_point(seg_desc.get_start_x() + seg_desc.get_length_x() * static_cast<double>(j) / (static_cast<double>(seg_desc.get_points_x()) - 1.0));

            // check segment_boundary points, since they might be in the mesh already:
            if (j == 0 || j == seg_desc.get_points_x() - 1)
            {
              // find matching segment boundary point:
              bool found = false;
              for (std::size_t k=0; k<segment_boundary_points.size(); ++k)
              {
                if (viennagrid::norm_2(candidate_point - segment_boundary_points[k]) < distance_tolerance )
                {
                  if (segment_boundary_ids[k] == -1) // point hasn't been added yet, so add now:
                  {
                    viennagrid::make_vertex_with_id(mesh, typename VertexType::id_type(vertex_counter), candidate_point);
                    segment_boundary_ids[k] = vertex_counter++;
                  }
                  segment_vertex_ids[i][j] = segment_boundary_ids[k];
                  found = true;
                  break;
                }
              }
              assert(found && bool("Logic error: Could not find matching boundary point"));
              (void)found; //silence unused variable warning in release mode
            }
            else
            {
              viennagrid::make_vertex_with_id(mesh, typename VertexType::id_type(vertex_counter), candidate_point);
              segment_vertex_ids[i][j] = vertex_counter++;
            }
          }
        }


        //
        // Set up cells:
        //
        log::info<log_generate_device>() << "* generate_device(): Setting up cells..." << std::endl;
        viennagrid::static_array<VertexHandleType, viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num> cell_vertex_handles;
        int cell_id = 0;
        for (std::size_t i=0; i<conf.size(); ++i)
        {
          SegmentDescriptionType const & seg_desc = conf.at(i);

          for (std::size_t j = 0; j<seg_desc.get_points_x() - 1; ++j)
          {
            cell_vertex_handles[0] = viennagrid::vertices(mesh).handle_at(static_cast<std::size_t>(segment_vertex_ids[i][j]));
            cell_vertex_handles[1] = viennagrid::vertices(mesh).handle_at(static_cast<std::size_t>(segment_vertex_ids[i][j+1]));

            viennagrid::make_element_with_id<CellType>(segmentation[static_cast<int>(i)],
                                                       cell_vertex_handles.begin(),
                                                       cell_vertex_handles.end(),
                                                       typename CellType::id_type(cell_id++));
          }
        }

        log::info<log_generate_device>() << "* generate_device(): DONE" << std::endl;
      }

      /** @brief Compatibility overload for simplex-lines (viennagrid::simplex_tag<1> and viennagrid::hypercube_tag<1> are the same)
       */
      template <typename MeshType>
      void generate_device_impl(MeshType & mesh,
                                device_generation_config const & conf,
                                viennagrid::hypercube_tag<1>
                               )
      {
        generate_device_impl(mesh, conf, viennagrid::simplex_tag<1>());
      }


      //
      // 2d generation
      //

      /** @brief Implementation of a simple two-dimensional 'mesh' generation for rectangular grids. Fills a ViennaGrid mesh.
       *
       * @param mesh           An empty ViennaGrid mesh
       * @param segmentation   A mesh segmentation for setting up the segments
       * @param conf           The generation configuration object
       */
      template <typename MeshT, typename SegmentationT>
      void generate_device_impl(MeshT & mesh,
                                SegmentationT & segmentation,
                                device_generation_config const & conf,
                                viennagrid::quadrilateral_tag
                               )
      {
        typedef typename viennagrid::result_of::point<MeshT>::type     PointType;
        typedef typename viennagrid::result_of::vertex<MeshT>::type    VertexType;
        typedef typename viennagrid::result_of::cell<MeshT>::type      CellType;

        typedef typename viennagrid::result_of::cell_tag<MeshT>::type  CellTag;

        typedef typename viennagrid::result_of::handle<MeshT, viennagrid::vertex_tag>::type  VertexHandleType;

        typedef typename device_generation_config::segment_description   SegmentDescriptionType;


        //
        // Prepare vertices at segment boundaries (note that this is O(N^2) with respect to the number of segments N. Not expected to hurt in practice, though...):
        //
        std::vector<PointType> segment_boundary_points;

        for (std::size_t seg_idx=0; seg_idx<conf.size(); ++seg_idx)
        {
          SegmentDescriptionType const & seg_desc = conf.at(seg_idx);

          assert(seg_desc.get_points_x() > 1   && bool("Logic error: Not enough points in x-direction provided for segment"));
          assert(seg_desc.get_points_y() > 1   && bool("Logic error: Not enough points in y-direction provided for segment"));
          assert(seg_desc.get_length_x() > 0.0 && bool("Logic error: x-coordinate is degenerate in device generation."));
          assert(seg_desc.get_length_y() > 0.0 && bool("Logic error: y-coordinate is degenerate in device generation."));

          double distance_tolerance = 1e-10 * std::min(seg_desc.get_length_x(), seg_desc.get_length_y());

          for (std::size_t j=0; j<seg_desc.get_points_y(); ++j)
          {
            for (std::size_t i=0; i<seg_desc.get_points_x(); ++i)
            {
              if (i == 0 || j == 0 || i == seg_desc.get_points_x() - 1 || j == seg_desc.get_points_y() - 1)
              {
                double inc_x = seg_desc.get_length_x() / static_cast<double>(seg_desc.get_points_x() - 1);
                double inc_y = seg_desc.get_length_y() / static_cast<double>(seg_desc.get_points_y() - 1);

                PointType p(seg_desc.get_start_x() + i * inc_x,
                            seg_desc.get_start_y() + j * inc_y);

                bool insert_p = true;

                for (std::size_t k=0; k<segment_boundary_points.size(); ++k)
                {
                  if (viennagrid::norm_2(p - segment_boundary_points[k]) < distance_tolerance)
                  {
                    insert_p = false;
                    break;
                  }
                }

                if (insert_p)
                  segment_boundary_points.push_back(p);
              }
            }
          }
        }

        //
        // Set up the vertices:
        //
        log::info<log_generate_device>() << "* generate_device(): Setting up vertices..." << std::endl;

        std::vector<long> segment_boundary_ids(segment_boundary_points.size(), -1);

        std::vector< std::vector< std::vector<long> > > segment_vertex_ids(conf.size()); // store the global vertex ID for each segment to save lookups later on.

        int vertex_counter = 0;

        // iterate over all segments and insert points if not already inserted in mesh. Store vertex ID for each segment.
        for (std::size_t seg_idx=0; seg_idx<conf.size(); ++seg_idx)
        {
          SegmentDescriptionType const & seg_desc = conf.at(seg_idx);

          segment_vertex_ids[seg_idx].resize(seg_desc.get_points_x());

          double distance_tolerance = 1e-10 * seg_desc.get_length_x();

          for (std::size_t i = 0; i<seg_desc.get_points_x(); ++i)
            segment_vertex_ids[seg_idx][i].resize(seg_desc.get_points_y());

          // Add get_points_x() points in the segment to the mesh if they haven't been added yet (check segment boundaries)
          double inc_x = seg_desc.get_length_x() / static_cast<double>(seg_desc.get_points_x() - 1);
          double inc_y = seg_desc.get_length_y() / static_cast<double>(seg_desc.get_points_y() - 1);

          for (std::size_t j = 0; j<seg_desc.get_points_y(); ++j)
          {
            for (std::size_t i = 0; i<seg_desc.get_points_x(); ++i)
            {
              PointType p(seg_desc.get_start_x() + i * inc_x,
                          seg_desc.get_start_y() + j * inc_y);

              // check segment_boundary points, since they might be in the mesh already:
              if (i == 0 || j == 0 || i == seg_desc.get_points_x() - 1 || j == seg_desc.get_points_y() - 1)
              {
                // find matching segment boundary point:
                bool found = false;
                for (std::size_t k=0; k<segment_boundary_points.size(); ++k)
                {
                  if (viennagrid::norm_2(p - segment_boundary_points[k]) < distance_tolerance )
                  {
                    if (segment_boundary_ids[k] == -1) // point hasn't been added yet, so add now:
                    {
                      viennagrid::make_vertex_with_id(mesh, typename VertexType::id_type(vertex_counter), p);
                      segment_boundary_ids[k] = vertex_counter++;
                    }
                    segment_vertex_ids[seg_idx][i][j] = segment_boundary_ids[k];
                    found = true;
                    break;
                  }
                }
                assert(found && bool("Logic error: Could not find matching boundary point"));
                (void)found; //silence unused variable warning in release mode
              }
              else
              {
                viennagrid::make_vertex_with_id(mesh, typename VertexType::id_type(vertex_counter), p);
                segment_vertex_ids[seg_idx][i][j] = vertex_counter++;
              }
            }
          }
        }


        //
        // Set up cells:
        //
        log::info<log_generate_device>() << "* generate_device(): Setting up cells..." << std::endl;
        viennagrid::static_array<VertexHandleType, viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num> cell_vertex_handles;
        int cell_id = 0;
        for (std::size_t seg_idx=0; seg_idx<conf.size(); ++seg_idx)
        {
          SegmentDescriptionType const & seg_desc = conf.at(seg_idx);

          for (std::size_t j = 0; j<seg_desc.get_points_y() - 1; ++j)
          {
            for (std::size_t i = 0; i<seg_desc.get_points_x() - 1; ++i)
            {
              cell_vertex_handles[0] = viennagrid::vertices(mesh).handle_at(std::size_t(segment_vertex_ids[seg_idx][i][j]));
              cell_vertex_handles[1] = viennagrid::vertices(mesh).handle_at(std::size_t(segment_vertex_ids[seg_idx][i+1][j]));
              cell_vertex_handles[2] = viennagrid::vertices(mesh).handle_at(std::size_t(segment_vertex_ids[seg_idx][i][j+1]));
              cell_vertex_handles[3] = viennagrid::vertices(mesh).handle_at(std::size_t(segment_vertex_ids[seg_idx][i+1][j+1]));

              viennagrid::make_element_with_id<CellType>(segmentation[static_cast<int>(seg_idx)],
                                                         cell_vertex_handles.begin(),
                                                         cell_vertex_handles.end(),
                                                         typename CellType::id_type(cell_id++));
            }
          }
        }

        log::info<log_generate_device>() << "* generate_device(): DONE" << std::endl;
      }


    } // namespace detail

    //
    // Public interface
    //

    /** @brief Public interface for mesh generation of simple one- or two-dimensional meshes. Device must be rectangular for the two-dimensional case.
     *
     * @param mesh    An empty ViennaGrid mesh
     * @param seg     An empty ViennaGrid segementation
     * @param conf    The mesh generation configuration object
     */
    template <typename MeshT, typename SegmentationT>
    void generate_device(MeshT & mesh, SegmentationT & seg, device_generation_config const & conf)
    {
      viennashe::util::detail::generate_device_impl(mesh, seg, conf, typename viennagrid::result_of::cell_tag<MeshT>::type());
    }

    /**
     * @brief A device generator to generate a device from C-Arrays (DOES NOT TAKE OWNERSHIP)
     * @param vertices A 2d array: vertices[num_vertices][dimension]
     * @param cells A 2d array: cells[num_cells][number_of_vertices_per_cell]
     * @param num_vertices The number of vertices in the domain
     * @param num_cells The number of cells in the domain
     */
    template < typename IndexT = unsigned long>
    struct device_from_array_generator
    {
      device_from_array_generator( double ** vertices,
                                   IndexT ** cells,
                                   IndexT * segmentation,
                                   IndexT   num_vertices,
                                   IndexT   num_cells)
      : vertices_(vertices), cells_(cells), segmentation_(segmentation),
        num_vertices_(num_vertices), num_cells_(num_cells)
      { }

      /** @brief Functor interface. The mesh and the segmentation need to be given. */
      template < typename MeshT, typename SegmentationT >
      void operator()(MeshT & mesh, SegmentationT & seg) const
      {
        typedef typename viennagrid::result_of::point<MeshT>::type    PointType;
        typedef typename viennagrid::result_of::vertex<MeshT>::type   VertexType;
        typedef typename viennagrid::result_of::cell<MeshT>::type     CellType;

        typedef typename viennagrid::result_of::cell_tag<MeshT>::type  CellTag;

        typedef typename viennagrid::result_of::handle<MeshT, viennagrid::vertex_tag>::type  VertexHandleType;

        for (int i = 0; i < static_cast<int>(num_vertices_); ++i)
        {
          PointType p;
          for (std::size_t j = 0; j < (std::size_t )PointType::dim; ++j)
            p[j] = vertices_[i][j];

          viennagrid::make_vertex_with_id( mesh, typename VertexType::id_type(i), p );
        }

        for (int i = 0; i < static_cast<int>(num_cells_); ++i )
        {
          viennagrid::static_array<VertexHandleType, viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num> cell_vertex_handles;
          for (std::size_t j=0; j < (std::size_t )viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num; ++j)
          {
            cell_vertex_handles[j] = viennagrid::vertices(mesh).handle_at(cells_[i][j]);
          }

          IndexT seg_id = segmentation_[i];
          if (segmentation_ != 0) seg_id = segmentation_[i];

          viennagrid::make_element_with_id<CellType>(seg[static_cast<int>(seg_id)],
                                                     cell_vertex_handles.begin(),
                                                     cell_vertex_handles.end(),
                                                     typename CellType::id_type(i));
        }
      }

    private:
      double ** vertices_;
      IndexT ** cells_;
      IndexT *  segmentation_;
      IndexT    num_vertices_;
      IndexT    num_cells_;

    };

    /**
     * @brief A device generator to generate a device from <b>flat</b> C-Arrays (DOES NOT TAKE OWNERSHIP)
     * @param vertices A 1d array: vertices[num_vertices * dimension]
     * @param cells A 1d array: cells[num_cells * number_of_vertices_per_cell]
     * @param num_vertices The number of vertices in the domain
     * @param num_cells The number of cells in the domain
     */
    template < typename IndexT = unsigned long>
    struct device_from_flat_array_generator
    {
      device_from_flat_array_generator( double * vertices,
                                        IndexT * cells,
                                        IndexT * segmentation,
                                        IndexT num_vertices,
                                        IndexT num_cells)
      : vertices_(vertices), cells_(cells), segmentation_(segmentation),
        num_vertices_(num_vertices), num_cells_(num_cells)
      { }

      /** @brief Functor interface. The mesh and the segmentation need to be given. */
      template < typename MeshT, typename SegmentationT >
      void operator()(MeshT & mesh, SegmentationT & seg) const
      {
        typedef typename viennagrid::result_of::point<MeshT>::type    PointType;
        typedef typename viennagrid::result_of::vertex<MeshT>::type   VertexType;
        typedef typename viennagrid::result_of::cell<MeshT>::type     CellType;

        typedef typename viennagrid::result_of::cell_tag<MeshT>::type  CellTag;

        typedef typename viennagrid::result_of::handle<MeshT, viennagrid::vertex_tag>::type  VertexHandleType;

        const std::size_t dim = (std::size_t )PointType::dim;
        for (int i = 0; i < static_cast<int>(num_vertices_); ++i)
        {
          PointType p;
          for (std::size_t j = 0; j < dim; ++j)
            p[j] = vertices_[std::size_t(i)*dim + j];

          viennagrid::make_vertex_with_id( mesh, typename VertexType::id_type(i), p );
        }

        for (int i = 0; i < static_cast<int>(num_cells_); ++i )
        {
          const std::size_t nen = viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num;

          viennagrid::static_array<VertexHandleType, viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num> cell_vertex_handles;
          for (std::size_t j = 0; j < nen; ++j)
          {
            cell_vertex_handles[j] = viennagrid::vertices(mesh).handle_at(cells_[std::size_t(i)*nen + j]);
          }

          IndexT seg_id = 0;
          if (segmentation_ != 0) seg_id = segmentation_[i];

          viennagrid::make_element_with_id<CellType>(seg[int(seg_id)],
                                                     cell_vertex_handles.begin(),
                                                     cell_vertex_handles.end(),
                                                     typename CellType::id_type(i));
        }
      }

    private:
      double        * vertices_;
      IndexT * cells_;
      IndexT * segmentation_;
      IndexT   num_vertices_;
      IndexT   num_cells_;

    };

  } //namespace util
} //namespace viennashe
#endif
