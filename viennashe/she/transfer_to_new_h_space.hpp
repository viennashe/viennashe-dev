#ifndef VIENNASHE_SHE_TRANSFER_TO_NEW_H_SPACE_HPP
#define VIENNASHE_SHE_TRANSFER_TO_NEW_H_SPACE_HPP

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

#include <cmath>

// viennagrid
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// viennashe
#include "viennashe/math/constants.hpp"
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/dispersion.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/she/harmonics_coupling.hpp"
#include "viennashe/she/assemble_common.hpp"
#include "viennashe/she/postproc/carrier_density.hpp"
#include "viennashe/she/postproc/current_density.hpp"

#include "viennashe/util/block_matrix_writer.hpp"
#include "viennashe/util/checks.hpp"
#include "viennashe/util/misc.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/she/timestep_quantities.hpp"

/** @file viennashe/she/transfer_to_new_h_space.hpp
    @brief Interpolation of a solution defined in one H-space to another H-space (required by Newton's method)
*/

namespace viennashe
{
  namespace she
  {

    namespace detail
    {

      template <typename ElementType,
                typename VertexT, typename EdgeT>
      void transfer_to_new_H_space_on_element(ElementType const & el,
                                              std::size_t index_H,
                                              she::unknown_she_quantity<VertexT, EdgeT> const & old_quan,
                                              she::unknown_she_quantity<VertexT, EdgeT> & new_quan)
      {
        if (new_quan.get_unknown_index(el, index_H) < 0) //nothing to do
          return;

        double new_kin_energy = new_quan.get_kinetic_energy(el, index_H);

        long old_index_H = energy_index_H(old_quan, el, static_cast<long>(std::min<std::size_t>(index_H, old_quan.get_value_H_size()-1)), new_kin_energy);

        if (old_index_H >= 0)  //old quantity does not have this kinetic energy, so set values to zero on new quantity
        {
          if (old_quan.get_expansion_order(el, std::size_t(old_index_H)) == 0)  //interpolate towards band-edge
          {
            while (old_quan.get_expansion_order(el, std::size_t(old_index_H)) == 0 && old_index_H < static_cast<long>(old_quan.get_value_H_size()) - 1) //box might be just below the band edge
              ++old_index_H;
          }

          if (old_quan.get_expansion_order(el, std::size_t(old_index_H)) > 0)  //there are values defined
          {
            std::size_t max_unknown_num = std::max<std::size_t>(new_quan.get_unknown_num(el,                 index_H ),
                                                                old_quan.get_unknown_num(el, std::size_t(old_index_H)));
            std::vector<double> values(max_unknown_num);
            for (std::size_t i=0; i<old_quan.get_unknown_num(el, std::size_t(old_index_H)); ++i)
              values[i] = old_quan.get_values(el, std::size_t(old_index_H))[i];
            new_quan.set_values(el, index_H, &(values[0]));
          }
        }

      }

      template <typename ElementType,
                typename VertexT, typename EdgeT>
      void normalize_on_new_H_space(ElementType const & el,
                                    she::unknown_she_quantity<VertexT, EdgeT> & quan,
                                    double scaling_factor)
      {
        for (std::size_t index_H=0; index_H < quan.get_value_H_size(); ++index_H)
        {
          if (quan.get_unknown_index(el, index_H) >= 0)
          {
            std::vector<double> values(quan.get_unknown_num(el, index_H));
            for (std::size_t i=0; i<values.size(); ++i)
              values[i] = quan.get_values(el, index_H)[i] * scaling_factor;
            quan.set_values(el, index_H, &(values[0]));
          }
        }
      }

    } // namespace detail


    /** @brief Interface transferring a solution given by 'old_solution' on some other (old) grid to the new grid.
    *
    * @param device         The device on which simulation is carried out
    * @param conf           The simulator configuration
    * @param new_quantities   The timestep_quantities used for the upcoming simulation, which is filled with the interpolated old values
    * @param old_quantities   The timestep_quantities used for obtaining the old values from another grid (same grid in (x, H)-space, but other potential profile)
    */
    template <typename DeviceType>
    void transfer_to_new_h_space(DeviceType const & device,
                                 viennashe::she::timestep_quantities<DeviceType> const & old_quantities,
                                 viennashe::she::timestep_quantities<DeviceType> & new_quantities,
                                 viennashe::config const & conf)
    {
      typedef typename DeviceType::mesh_type              MeshType;

      typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;
      typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;

      typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
      typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;

      typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

      typedef she::unknown_she_quantity<CellType, FacetType>     SHEQuantity;

      MeshType const & mesh = device.mesh();

      for (std::size_t i=0; i<new_quantities.unknown_she_quantities().size(); ++i)
      {
        SHEQuantity const & old_quan = old_quantities.unknown_she_quantities()[i];
        SHEQuantity       & new_quan = new_quantities.unknown_she_quantities()[i];

        if (new_quan.get_value_H_size() == 0)  //don't do anything if quantity is not enabled
          return;

        viennashe::she::carrier_density_wrapper<SHEQuantity> old_density_wrapper(conf, old_quan);
        viennashe::she::carrier_density_wrapper<SHEQuantity> new_density_wrapper(conf, new_quan);

        viennashe::she::detail::current_on_facet_by_ref_calculator<DeviceType, SHEQuantity>
           old_edge_evaluator(device, conf, old_quan);
        viennashe::she::detail::current_on_facet_by_ref_calculator<DeviceType, SHEQuantity>
           new_edge_evaluator(device, conf, new_quan);

        //
        // Step 1: transfer solution on vertices
        //
        CellContainer cells(mesh);
        for (CellIterator cit = cells.begin();
            cit != cells.end();
            ++cit)
        {
          for (std::size_t index_H = 1; index_H < new_quan.get_value_H_size() - 1; ++index_H)
            detail::transfer_to_new_H_space_on_element(*cit, index_H, old_quan, new_quan);

          // scale with respect to density:
          double density_old = old_density_wrapper(*cit);
          double density_new = new_density_wrapper(*cit);

          if (density_new > 0 && density_old > 0)
            detail::normalize_on_new_H_space(*cit, new_quan, density_old / density_new);
        } //for vertices


        //
        // Step 2: transfer solution on edges
        //
        FacetContainer facets(mesh);
        for (FacetIterator fit = facets.begin();
             fit != facets.end();
             ++fit)
        {
          for (std::size_t index_H = 1; index_H < new_quan.get_value_H_size()-1; ++index_H)
            detail::transfer_to_new_H_space_on_element(*fit, index_H, old_quan, new_quan);

          // scale with respect to current:
          double current_old = old_edge_evaluator(*fit);
          double current_new = new_edge_evaluator(*fit);

          if (current_new && current_old)
            detail::normalize_on_new_H_space(*fit, new_quan, current_old / current_new);
        }
      }


    } //transfer_to_new_H_space

  } //namespace she
} //namespace viennashe

#endif

