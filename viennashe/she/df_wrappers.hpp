#ifndef VIENNASHE_SHE_DF_WRAPPERS_HPP
#define VIENNASHE_SHE_DF_WRAPPERS_HPP

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

#include "viennashe/forwards.h"
#include "viennashe/physics/constants.hpp"
#include "viennashe/math/spherical_harmonics.hpp"

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/mesh.hpp"

#include "viennashe/config.hpp"
#include "viennashe/she/exception.hpp"
#include "viennashe/she/she_quantity.hpp"
#include "viennashe/simulator_quantity.hpp"
#include "viennashe/util/filter.hpp"
#include "viennashe/util/misc.hpp"

/** @file viennashe/she/df_wrappers.hpp
    @brief Provides a convenience wrapper for accessing the distribution function coefficients over the device.
*/

namespace viennashe
{
  namespace she
  {
    namespace detail
    {
      /** @brief Returns the best total energy in order to match the supplied kinetic energy */
      template <typename UnknownSHEType, typename ElementType>
      std::size_t find_best_H(UnknownSHEType const & she_quantity, ElementType const & el, double kin_energy, std::size_t index_H_guess)
      {
        long carrier_sign = (she_quantity.get_carrier_type_id() == viennashe::ELECTRON_TYPE_ID) ? -1 : 1;

        std::size_t index_H = std::min<std::size_t>(index_H_guess, she_quantity.get_value_H_size() - 2);

        // find best index_H:
        while (   (she_quantity.get_kinetic_energy(el, index_H) < kin_energy)  //search upwards
                && ( (carrier_sign == 1) ? (index_H > 0) : (index_H < she_quantity.get_value_H_size() - 1) )
              )
        {
          index_H -= static_cast<std::size_t>(carrier_sign);
        }

        while (   (she_quantity.get_kinetic_energy(el, index_H) > kin_energy)  //search downwards
                && ( (carrier_sign == -1) ? (index_H > 0) : (index_H < she_quantity.get_value_H_size() - 1) )
              )
        {
          index_H += static_cast<std::size_t>(carrier_sign);
        }

        // finally, pick index_H or index_H + 1 (electrons), whichever is better
        double error_down = she_quantity.get_kinetic_energy(el, index_H) - kin_energy;
        if (error_down <= 0)
          return index_H;

        double error_up   = error_down;
        if (   static_cast<long>(index_H) >= static_cast<long>(carrier_sign)
            && static_cast<long>(index_H) <  static_cast<long>(she_quantity.get_value_H_size()) + static_cast<long>(carrier_sign))
        {
          error_up = she_quantity.get_kinetic_energy(el, static_cast<std::size_t>(static_cast<long>(index_H) - carrier_sign)) - kin_energy;
          if (she_quantity.get_unknown_index(el, static_cast<std::size_t>(static_cast<long>(index_H) - carrier_sign)) < 0)
            error_up = 2.0 * error_down; //make sure error_up is not selected
        }

        if (std::abs(error_up) < std::abs(error_down))
          index_H -= static_cast<std::size_t>(carrier_sign);

        return index_H;
      }
    }



    //  Class                                                Arguments                         Comments
    // --------------------------------------------------------------------------------------------------------------------------------------------
    // she_df_wrapper                                        (tag, energy, l, m)               Returns f_{l,m}. Even l on vertices, odd l on edges. Optional guess argument
    //  |
    //  |--> interpolated_she_df_wrapper                     (tag, energy, l, m)               Returns f_{l,m} on vertices. Odd l are automatically interpolated from edges.
    //  |      |
    //         \--> df_wrapper  --> generalized_df_wrapper   (tag, energy, theta, phi)         Returns the full (generalized) distribution function on each vertex.
    //  |
    //  \---------> edf_wrapper --> generalized_edf_wrapper  (tag, energy)                     Returns the (generalized) energy distribution function on each vertex.
    //

    /** @brief A convenience wrapper for accessing the distribution function coefficients over the device. Even-order coefficients can be obtained from vertices, odd-order coefficients from edges.
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class she_df_wrapper
    {
      private:
        typedef she_df_wrapper self_type;

        typedef typename DeviceType::mesh_type           MeshType;

        typedef typename viennagrid::result_of::vertex<MeshType>::type       VertexType;
        typedef typename viennagrid::result_of::facet<MeshType>::type        FacetType;
        typedef typename viennagrid::result_of::cell<MeshType>::type         CellType;

      public:

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;


        /** @brief Constructs the wrapper.
         *
         * @param quan             The SHE quantities
         * @param conf             The simulator configuration
        */
        she_df_wrapper(viennashe::config const & conf,
                       SHEQuantityT const & quan)
          : conf_(conf), she_unknown_(quan), dispersion_(conf.dispersion_relation(quan.get_carrier_type_id()))
        { }


        /** @brief Returns an even-order SHE coefficient for the respective carrier on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param l               Major (leading) spherical harmonics index
         *  @param m               Minor spherical harmonics index
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          std::size_t l,
                          long m,
                          std::size_t index_H_guess = 0) const
        {
          if ( (l % 2 != 0) || (std::abs(static_cast<double>(m)) > static_cast<double>(l)) )
            throw invalid_expansion_order_exception("Invalid expansion order on vertex", l, m);

          return evaluate(cell,
                          detail::find_best_H(she_unknown_, cell, kinetic_energy, index_H_guess),
                          l, m);
        }


        /** @brief Returns an even-order SHE coefficient for the respective carrier on an edge:
         *
         *  @param facet           The facet from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param l               Major (leading) spherical harmonics index
         *  @param m               Minor spherical harmonics index
         *  @param index_H_guess   Initial guess for the lookup of index_H (optimization purposes...)
         */
        double operator()(FacetType const & facet,
                          double kinetic_energy,
                          std::size_t l,
                          long m,
                          std::size_t index_H_guess = 0) const
        {
          if ( (l % 2 != 1) || (std::abs(static_cast<double>(m)) > static_cast<double>(l)) )
            throw invalid_expansion_order_exception("Invalid expansion order on edge", l, m);

          return evaluate(facet,
                          detail::find_best_H(she_unknown_, facet, kinetic_energy, index_H_guess),
                          l, m);
        }


        /** @brief Batch-evaluation: Returns the expansion coefficents computed at the particular vertex for the particular kinetic energy */
        template <typename FoldVectorType>
        void fill(CellType const & cell,
                  double kin_energy,
                  std::size_t index_H_guess,
                  FoldVectorType & vec)  const
        {
          typedef typename FoldVectorType::value_type   value_type;

          std::size_t index_H = detail::find_best_H(she_unknown_, cell, kin_energy, index_H_guess);

          value_type const * pValues = she_unknown_.get_values(cell, index_H);

          std::size_t num_values = she_unknown_.get_unknown_num(cell, index_H);
          for (std::size_t i=0; i < num_values; ++i)
            vec[i] = pValues[i];
        }


        /** @brief Batch-evaluation: Returns the expansion coefficents computed at the particular edge for the particular kinetic energy */
        template <typename FoldVectorType>
        void fill(FacetType const & facet,
                  double kin_energy,
                  std::size_t index_H_guess,
                  FoldVectorType & vec) const
        {
          typedef typename FoldVectorType::value_type   value_type;

          std::size_t index_H = detail::find_best_H(she_unknown_, facet, kin_energy, index_H_guess);

          value_type const * pValues = she_unknown_.get_values(facet, index_H);

          std::size_t num_values = she_unknown_.get_unknown_num(facet, index_H);
          for (std::size_t i=0; i < num_values; ++i)
            vec[i] = pValues[i];
        }

        viennashe::config const & config() const { return conf_; }
        SHEQuantityT const & quan() const { return this->she_unknown_; }

        dispersion_relation_type const & dispersion_relation() const { return this->dispersion_; }

      private:

        /** @brief Implementation of the distribution function evaluation
         *
         * @param el         Either a vertex or an edge
         * @param index_H    Total energy index
         * @param l          Leading spherical harmonics index
         * @param m          Minor spherical harmonics index
         * @param parity     0: even order, 1: odd order
         */
        template <typename ElementType>
        double evaluate(ElementType const & el,
                        std::size_t index_H,
                        std::size_t l,
                        long m) const
        {
          std::size_t vec_index = 0;
          const std::size_t unknown_num = she_unknown_.get_unknown_num(el, index_H);

          for (std::size_t i = parity_for_element(el); i < l; i += 2)
            vec_index += 2 * i + 1;

          if (unknown_num > 0 && unknown_num >= vec_index)
          {
            switch (conf_.she_discretization_type())
            {
            case SHE_DISCRETIZATION_EVEN_ODD_ORDER_DF:
              return she_unknown_.get_values(el, index_H)[vec_index + std::size_t(long(l) + m)];
            case SHE_DISCRETIZATION_EVEN_ODD_ORDER_GENERALIZED_DF:
              return she_unknown_.get_values(el, index_H)[vec_index + std::size_t(long(l) + m)] / averaged_density_of_states(she_unknown_, dispersion_, el, index_H);
            default: throw std::runtime_error("she_df_wrapper::evaluate(): Unknown SHE discretization type!");
            }
          }

          return 0.0;
        }

        std::size_t parity_for_element(CellType  const &) const { return 0; }
        std::size_t parity_for_element(FacetType const &) const { return 1; }

        viennashe::config conf_;
        SHEQuantityT she_unknown_;
        dispersion_relation_type  dispersion_;
    };


    //
    // Start of derived wrappers: Interpolated SHE wrapper branch
    //

    /** @brief A convenience wrapper for accessing the distribution function coefficients over the device on each vertex. Note that odd-order expansion coefficients are averaged, thus one MUST NOT evaluate flux-quantities using this interpolated wrapper!
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class interpolated_she_df_wrapper
    {
      private:
        typedef typename DeviceType::mesh_type           MeshType;

        typedef typename viennagrid::result_of::facet<MeshType>::type      FacetType;
        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      public:

        typedef interpolated_she_df_wrapper self_type;

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;

        interpolated_she_df_wrapper(DeviceType const & device,
                                    viennashe::config const & conf,
                                    SHEQuantityT const & quan)
          : device_(device),      //Note: Reference to device is required for ViennaGrid co-boundary operations
            she_df_(conf, quan) {}

        /** @brief Returns the energy distribution function for the respective carrier type in one valley on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param l               Major (leading) spherical harmonics index
         *  @param m               Minor spherical harmonics index
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          std::size_t l,
                          long m,
                          std::size_t index_H_guess = 0) const
        {
          if (l % 2 == 0)
            return she_df_(cell, kinetic_energy, l, m, index_H_guess);

          return interpolated_odd_coefficient(cell, kinetic_energy, l, m, index_H_guess);
        }


        SHEQuantityT const & quan() const { return this->she_df_.quan(); }

        dispersion_relation_type const & dispersion_relation() const { return this->she_df_.dispersion_relation(); }

      private:

        double interpolated_odd_coefficient(CellType const & cell,
                                            double kinetic_energy,
                                            std::size_t l,
                                            long m,
                                            std::size_t index_H_guess) const
        {
          typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;
          typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type  FacetOnCellIterator;

          // Interpolation strategy: Due to the Galerkin approach for the odd-order coefficients, the weighted arithmetic mean is returned. Dual box interpolator does not make sense here, because this is an interpolation from scalar to scalar.
          double result = 0;
          double summed_volume = 0;

          FacetOnCellContainer facets(cell);
          for (FacetOnCellIterator focit = facets.begin();
                                   focit != facets.end();
                                 ++focit)
          {
            // TODO: Think about how to deal with facets which don't carry unknowns!

            double facet_volume = viennagrid::volume(*focit);
            summed_volume += facet_volume;
            result += she_df_(*focit, kinetic_energy, l, m, index_H_guess) * facet_volume;
            //Note: Alternatives for interpolations are:
            //      - Harmonic weights: 'result += she_df_() / box_volume; ...; result *= summed_volume;'
            //      - No weights: Just ignore box_volume.
          }

          if (summed_volume > 0)
            result /= summed_volume;

          return result;
        }

        DeviceType const & device_;
        she_df_wrapper<DeviceType, SHEQuantityT>  she_df_;
    };



    /** @brief A convenience wrapper for evaluating the full distribution function.
     * Note: Evaluations are quite costly. If you wish to evaluate integrals over the distribution function,
     * you should definitely consider the use of SHE expansion coefficients as well as the orthogonality relations for spherical harmonics.
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class df_wrapper
    {
      private:
        typedef df_wrapper self_type;

        typedef typename DeviceType::mesh_type           MeshType;

        typedef typename viennagrid::result_of::facet<MeshType>::type      FacetType;
        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

        typedef unknown_she_quantity<CellType, FacetType>   UnknownSHEType;

      public:

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;

        df_wrapper(DeviceType const & device,
                   viennashe::config const & conf,
                   SHEQuantityT const & quan)
        : she_unknown_(quan),
          interpolated_she_df_(device, conf, quan) {}

        /** @brief Returns the energy distribution function for electrons in one valley on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param theta           Polar angle
         *  @param phi             Azimuthal angle
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          double theta,
                          double phi,
                          std::size_t index_H_guess = 0) const
        {
          return evaluate(cell, kinetic_energy, theta, phi, index_H_guess);
        }

        SHEQuantityT const & quan() const { return this->interpolated_she_df_.quan(); }

        dispersion_relation_type const & dispersion_relation() const { return this->interpolated_she_df_.dispersion_relation(); }

      private:

        double evaluate(CellType const & cell,
                        double kinetic_energy,
                        double theta,
                        double phi,
                        std::size_t index_H_guess) const
        {
          std::size_t index_H = detail::find_best_H(she_unknown_, cell, kinetic_energy, index_H_guess);
          std::size_t L = she_unknown_.get_expansion_order(cell, index_H);

          double result = 0;
          for (std::size_t l=0; l <= L; ++l)
          {
            for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
            {
              viennashe::math::SphericalHarmonic Y_lm(static_cast<int>(l),m);
              result += interpolated_she_df_(cell, kinetic_energy, l, m, index_H) * Y_lm(theta, phi);
            }
          }

          return result;
        }

        UnknownSHEType she_unknown_;
        interpolated_she_df_wrapper<DeviceType, SHEQuantityT>  interpolated_she_df_;
    };




    /** @brief A convenience wrapper for evaluating the generalized distribution function (i.e. f * Z, with density of states Z). Note: Evaluations are quite costly. If you wish to evaluate integrals over the distribution function, you should definitely consider the use of SHE expansion coefficients as well as the orthogonality relations for spherical harmonics.
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class generalized_df_wrapper
    {
      private:
        typedef generalized_df_wrapper                 self_type;

        typedef typename DeviceType::mesh_type       MeshType;

        typedef typename viennagrid::result_of::facet<MeshType>::type      FacetType;
        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      public:

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;

        generalized_df_wrapper(DeviceType const & device,
                               viennashe::config const & conf,
                               SHEQuantityT const & quan) : dispersion_(conf.dispersion_relation(quan.get_carrier_type_id())), df_(device, conf, quan) {}

        /** @brief Returns the energy distribution function for the respective carrier in one valley on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param theta           Polar angle
         *  @param phi             Azimuthal angle
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          double theta,
                          double phi,
                          std::size_t index_H_guess = 0) const
        {
          return df_(cell, kinetic_energy, theta, phi, index_H_guess) * dispersion_.density_of_states(kinetic_energy, theta, phi);
        }

        SHEQuantityT const & quan() const { return this->df_.quan(); }

        dispersion_relation_type const & dispersion_relation() const { return this->df_.dispersion_relation(); }


      private:
        dispersion_relation_type      dispersion_;
        df_wrapper<DeviceType, SHEQuantityT>   df_;
    };






    //
    // Derived wrappers: EDF Wrappers
    //

    /** @brief A convenience wrapper for accessing the energy distribution function on each vertex of the device.
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class edf_wrapper
    {
      private:
        typedef edf_wrapper                           self_type;

        typedef typename DeviceType::mesh_type      MeshType;

        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      public:

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;


        /** @brief Constructs the wrapper.
         *
         * @param quan             The SHE quantities.
         * @param conf             The simulator configuration
         */
        edf_wrapper(viennashe::config const & conf,
                    SHEQuantityT const & quan) : she_df_(conf, quan), Y_00_(0, 0) {}


        /** @brief Returns the energy distribution function for the respective carrier in one valley on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          std::size_t index_H_guess = 0) const
        {
          return she_df_(cell, kinetic_energy, 0, 0, index_H_guess) * Y_00_(0, 0);
        }

        SHEQuantityT const & quan() const { return this->she_df_.quan(); }

        dispersion_relation_type const & dispersion_relation() const { return this->she_df_.dispersion_relation(); }


      private:
        she_df_wrapper<DeviceType, SHEQuantityT>  she_df_;
        viennashe::math::SphericalHarmonic         Y_00_;
    };


    /** @brief A convenience wrapper for accessing the generalized energy distribution function (f * Z, with density of states Z) on each vertex of the device.
     *
     * @tparam DeviceType      The device type on which to evaluate the distribution function
     * @tparam VectorType      Vector type used for the SHE result vector.
     */
    template <typename DeviceType, typename SHEQuantityT>
    class generalized_edf_wrapper
    {
      private:
        typedef generalized_edf_wrapper                    self_type;

        typedef typename DeviceType::mesh_type           MeshType;

        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      public:

        typedef typename viennashe::config::dispersion_relation_type      dispersion_relation_type;
        typedef SHEQuantityT she_quantity_type;

        /** @brief Constructs the wrapper.
         *
         * @param quan             The SHE quantities
         * @param conf             The simulator configuration
         */
        generalized_edf_wrapper(viennashe::config const & conf,
                                SHEQuantityT const & quan)
          : dispersion_(conf.dispersion_relation(quan.get_carrier_type_id())), edf_(conf, quan) {}


        /** @brief Returns the energy distribution function for electrons in one valley on a vertex:
         *
         *  @param cell            The cell from which the SHE coefficient should be returned
         *  @param kinetic_energy  Kinetic energy
         *  @param index_H_guess   The total energy index for starting the search for the best total energy index available. Allows to speed-up evaluation considerably.
         */
        double operator()(CellType const & cell,
                          double kinetic_energy,
                          std::size_t index_H_guess = 0) const
        {
          return edf_(cell, kinetic_energy, index_H_guess) * dispersion_.density_of_states(kinetic_energy);
        }

        SHEQuantityT const & quan() const { return this->edf_.quan(); }

        dispersion_relation_type const & dispersion_relation() const { return this->edf_.dispersion_relation(); }


      private:
        dispersion_relation_type      dispersion_;
        edf_wrapper<DeviceType, SHEQuantityT>  edf_;
    };


  }
}

#endif
