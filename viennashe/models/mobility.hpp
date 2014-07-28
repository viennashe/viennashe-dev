#ifndef VIENNASHE_MODELS_MOBILITY_HPP
#define VIENNASHE_MODELS_MOBILITY_HPP

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

// viennashe
#include "viennashe/forwards.h"

#include "viennashe/models/mobility_parameters.hpp"
#include "viennashe/models/mobility_model.hpp"

namespace viennashe
{
  namespace models
  {
    /** @brief Compiletime evaluation namespace */
    namespace result_of
    {
      /** @brief Compiletime mobility-type getter */
      template < typename DeviceType >
      struct mobility_type
      {
        typedef typename viennashe::models::dd::mobility<DeviceType> type;
      };

    } // result_of

    /**
     * @brief Creates a new mobility model using the given parameters
     * @param device The device
     * @param params The mobility model parameters
     * @return A viennashe::models::dd::mobility model
     */
    template < typename DeviceType >
    viennashe::models::dd::mobility<DeviceType>
      create_mobility_model(DeviceType const & device,
                            viennashe::models::dd::mobility_paramters const & params)
    {
      return viennashe::models::dd::mobility<DeviceType>(device, params);
    }

    /**
     * @brief Returns a mobility model, which always yields the same mobility
     * @param device
     * @param mu The mobility you like the model to yield
     * @return A viennashe::models::dd::mobility model, where all submodels have been disabled
     */
    template < typename DeviceType >
    viennashe::models::dd::mobility<DeviceType> create_constant_mobility_model(DeviceType const & device, double mu)
    {
      viennashe::models::dd::mobility_paramters params;
      params.mu0 = mu;
      params.field.enabled    = false;
      params.lattice.enabled  = false;
      params.impurity.enabled = false;
      params.surface.enabled  = false;

      return viennashe::models::dd::mobility<DeviceType>(device, params);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class mobility_base
    {
      //typedef typename DeviceT::mesh_type    MeshType;

    public:
      //typedef typename viennagrid::result_of::point<MeshType>::type    PointType;
      //typedef typename viennagrid::result_of::facet<MeshType>::type    FacetType;
      //typedef typename viennagrid::result_of::cell<MeshType>::type     CellType;

      //virtual double evaluate(CellType  const & cell_0,
      //                        FacetType const & facet,
      //                        CellType  const & cell_1) const = 0;

      // TODO: Fix this super-ugly workaround
      virtual double evaluate(void const *cell_0,
                              void const *facet,
                              void const *cell_1) const = 0;

      virtual mobility_base * clone() const = 0;

      virtual ~mobility_base() {}
    };



    //template<typename DeviceT>
    class constant_mobility : public mobility_base
    {
      //typedef typename DeviceT::mesh_type    MeshType;

    public:
      //typedef typename viennagrid::result_of::point<MeshType>::type    PointType;
      //typedef typename viennagrid::result_of::facet<MeshType>::type    FacetType;
      //typedef typename viennagrid::result_of::cell<MeshType>::type     CellType;

      /** @brief Initialize the constant mobility model with the given mobility in SI units: m^2/(Vs) */
      constant_mobility(double val) : mu_(val) {}

      double evaluate(void const *cell_0,
                      void const *facet,
                      void const *cell_1) const
      {
        (void)cell_0; (void)facet; (void)cell_1;
        return mu_;
      }

      template<typename CellT, typename FacetT>
      double evaluate(CellT  const & cell_0,
                      FacetT const & facet,
                      CellT  const & cell_1) const
      {
        (void)cell_0; (void)facet; (void)cell_1;
        return mu_;
      }

      mobility_base * clone() const { return new constant_mobility(mu_); }

    private:
      double mu_;
    };


    template<typename DeviceT>
    class hcd_mobility : public mobility_base
    {
      typedef typename DeviceT::mesh_type    MeshType;

    public:
      typedef typename viennagrid::result_of::point<MeshType>::type    PointType;
      typedef typename viennagrid::result_of::facet<MeshType>::type    FacetType;
      typedef typename viennagrid::result_of::cell<MeshType>::type     CellType;

      /** @brief Initialize the constant mobility model with the given mobility in SI units: m^2/(Vs) */
      hcd_mobility(std::vector<double> const & nit_values, double val) : nit_values_(nit_values), mu_(val) {}

      double evaluate(void const *cell_0,
                      void const *facet,
                      void const *cell_1) const
      {
        return evaluate(*reinterpret_cast<CellType const *>(cell_0),
                        *reinterpret_cast<FacetType const *>(facet),
                        *reinterpret_cast<CellType const *>(cell_1));
      }

      double evaluate(CellType  const & cell_0,
                      FacetType const & facet,
                      CellType  const & cell_1) const
      {
        PointType facet_centroid = viennagrid::centroid(facet);
        double nit_density = (  nit_values_[static_cast<std::size_t>(cell_0.id().get())]
                              + nit_values_[static_cast<std::size_t>(cell_1.id().get())] ) / 2.0;

        facet_centroid[1] = 0.0;
        return mu_ / (1.0 + 1e-18 * nit_density * std::exp(std::abs(facet_centroid[1])/1e-8));
      }

      mobility_base * clone() const { return new hcd_mobility(nit_values_, mu_); }

    private:
      std::vector<double> const & nit_values_;
      double mu_;
    };


    /** @brief A proxy object for dealing with mobility implementations. Does NOT take ownership of the provided pointer! */
    class mobility_proxy
    {
    public:
      mobility_proxy(mobility_base const * ptr) : ptr_(ptr) {}

      /** @brief Returns the mobility evaluated at the facet between two cells.
       */
      template<typename CellT, typename FacetT>
      double evaluate(CellT  const & cell_0,
                      FacetT const & facet,
                      CellT  const & cell_1) const
      {
        return ptr_->evaluate(reinterpret_cast<void const *>(&cell_0),
                              reinterpret_cast<void const *>(&facet),
                              reinterpret_cast<void const *>(&cell_1));
      }

      mobility_base const * get() const { return ptr_; }

    private:
      mobility_base const * ptr_;
    };

  } // models
} // viennashe

#endif /* VIENNASHE_MODELS_MOBILITY_HPP */

