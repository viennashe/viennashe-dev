#ifndef VIENNASHE_FORWARDS_H
#define VIENNASHE_FORWARDS_H

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

/**
 @mainpage Main Page

 On these pages you can find all the documentation about the free open source semiconductor device simulator ViennaSHE.
 We recommend starting with the \ref manual-page-introduction, which provides a general overview of the primary solution method employed.

 ViennaSHE is designed in a library-centric way, building on top of:
   - the GPU-accelerated linear solver library [ViennaCL](http://viennacl.sourceforge.net/)
   - the generic mesh datastructure library [ViennaGrid](http://viennagrid.sourceforge.net/).

 A good starting point for your own simulations is to pick one of our [Examples](./examples.html) and alter it accordingly to fit your needs.

<i>User contributions welcome! If you wish to contribute a feature to ViennaSHE, do not hesitate to contact us e.g. via the mailing list at viennashe-support@lists.sourceforge.net</i>


 -----------------------------------
 \htmlonly
 <div style="align: right; width: 100%">
 <a href="http://www.tuwien.ac.at/"><img src="tuwien.png"></a>
 <a href="http://www.iue.tuwien.ac.at/"><img src="iue.png"></a>
 <a href="http://www.asc.tuwien.ac.at/"><img src="asc.png"></a>
 </div>
 \endhtmlonly
*/

#include <stdexcept>

#include <ostream>
#include <cstddef>
#include <string>
#include <map>

#include "viennashe/version.hpp"

#include "viennagrid/forwards.hpp"


/** @file viennashe/forwards.h
    @brief Contains forward declarations and definition of small classes that must be defined at an early stage
*/

//
//  Section 1: Set up keys for ViennaData:
//

/** @brief The main ViennaSHE namespace. All functionality resides inside this namespace. */
namespace viennashe
{
  /** @brief Namespace for implementation details. Typically not of interest for the library user. */
  namespace detail
  {
    /** @brief Metafunction for the determination of whether a mesh is 1d. */
    template <typename MeshT, typename CellTag = typename viennagrid::result_of::cell_tag<MeshT>::type>
    struct is_1d_mesh
    {
      enum { value = false };
    };

    /** @brief Metafunction for the determination of whether a mesh is 1d. Specialization for 1d. */
    template <typename MeshT>
    struct is_1d_mesh<MeshT, viennagrid::simplex_tag<1> >
    {
      enum { value = true };
    };

    /** @brief Metafunction for the determination of whether a mesh is 1d. Specialization for 1d. */
    template <typename MeshT>
    struct is_1d_mesh<MeshT, viennagrid::hypercube_tag<1> >
    {
      enum { value = true };
    };
  }

  /** @brief A helper class to raise compile time errors */
  template <typename T>
  struct error_indicator {};


  template <typename MeshT,
            bool edges_and_cells_different = detail::is_1d_mesh<MeshT>::value >
  class device;

  template <typename DeviceType>
  class simulator;


  /** @brief Information on the number of unknowns per quantity.
    *
    * First member of pair denotes the number of unknowns defined on vertices
    * Second member of pair denotes the number of unknowns defined on edges
    */
  typedef std::map<std::string, std::pair<std::size_t, std::size_t> >    map_info_type;


  /** @brief An enumeration of all equation types ViennaSHE supports */
  enum equation_id
  {
    EQUATION_INVALID = 0,
    EQUATION_POISSON_DD,
    EQUATION_CONTINUITY,
    EQUATION_SHE,
    EQUATION_DENSITY_GRADIENT,
    EQUATION_POISSON_HEAT
  };

  /** @brief An enumeration of all boundary conditions ViennaSHE supports */
  enum boundary_type_id
  {
    BOUNDARY_NONE = 0,
    BOUNDARY_DIRICHLET,      //     u = beta
    BOUNDARY_NEUMANN,        // du/dn = beta
    BOUNDARY_ROBIN,           // du/dn = beta - alpha * u
    BOUNDARY_GENERATION_RECOMBINATION // for SHE
  };

  enum material_category_id
  {
    MATERIAL_CONDUCTOR_ID,
    MATERIAL_NO_CONDUCTOR_ID,
    MATERIAL_SEMICONDUCTOR_ID,
    MATERIAL_NO_SEMICONDUCTOR_ID,
    MATERIAL_INSULATOR_ID,
    MATERIAL_NO_INSULATOR_ID
  };

  enum she_discretization_type_id
  {
    SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF,
    SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF
  };

  enum she_scaling_type_id
  {
    SHE_SCALING_KINETIC_ENERGY,
    SHE_SCALING_TOTAL_ENERGY
  };

  namespace math
  {
    /** @brief An enumeration of spherical harmonics types (all, even, odd) ViennaSHE supports */
    enum harmonics_iteration_type
    {
      ALL_HARMONICS_ITERATION_ID = 0,
      EVEN_HARMONICS_ITERATION_ID,
      ODD_HARMONICS_ITERATION_ID
    };
  }

  //
  // Box-integration related. Keys used for storing box integration information on the mesh, cf. ViennaGrid manual.
  //
  /** @brief Internal key for storing the box volume associated with an edge or a vertex */
  struct box_volume_key {};            //box volume associated with an edge or vertex

  /** @brief Internal key for storing the interface area of a box associated with an edge */
  struct edge_interface_area_key {};   //box volume associated with an edge


  // Other ViennaSHE forward declarations:
  /** @brief Type for storing the discrete energies */
  typedef std::vector<double>                   energy_vector_type;

  /** @brief Type for storing the unknown indices at each point in the (x, H)-space. */
  typedef std::vector<long>                     she_index_vector_type;

  /** @brief Enumeration type for selecting the carrier type */
  enum carrier_type_id
  {
    INVALID_TYPE = 0,
    ELECTRON_TYPE_ID,
    HOLE_TYPE_ID
  };

  /** @brief Enumeration type for selecting the doping type */
  enum doping_type_id
  {
    INVALID_DOPING_TYPE = 0,
    DONATOR_DOPING_TYPE_ID,
    ACCEPTOR_DOPING_TYPE_ID
  };

  // Solver tags:
  namespace solvers
  {
    /** @brief Internal tag used for the specification of a dense linear solver (Gauss, single-threaded) */
    class dense_linear_solver_tag {};

    /** @brief Internal tag used for the specification of a single-threaded linear solver */
    class serial_linear_solver_tag {};

    /** @brief Internal tag used for the specification of a CPU-based multi-threaded linear solver */
    class parallel_linear_solver_tag {};

    /** @brief Internal tag used for the specification of a CPU-based PETSC solver */
    class petsc_linear_solver_tag {};

    class petsc_amgx_solver_tag {};
    /** @brief Internal tag used for the specification of a GPU-accelerated linear solver */
    class gpu_parallel_linear_solver_tag {};

  }
}


#endif
