#ifndef VIENNASHE_SHE_SCATTERING_COMMON_HPP
#define VIENNASHE_SHE_SCATTERING_COMMON_HPP
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
#include "viennashe/math/constants.hpp"
#include "viennashe/materials/all.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/physics/physics.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"
#include "viennashe/config.hpp"

/** @file viennashe/she/scattering/common.hpp
    @brief Common classes for scattering operators
*/

namespace viennashe
{
  namespace she
  {

    enum scatter_process_id
    {
      INVALID_SCATTER_PROCESS = 0,
      ACOUSTIC_PHONON_SCATTERING,
      FIXED_CHARGE_SCATTERING,
      IMPACT_IONIZATION_SCATTERING,
      IMPURITY_SCATTERING,
      OPTICAL_PHONON_SCATTERING,
      SURFACE_ACOUSTIC_PHONON_SCATTERING,
      SURFACE_ROUGHNESS_SCATTERING,
      SURFACE_SCATTERING,
      TRAPPED_CHARGE_SCATTERING
    };


    /** @brief A simple class returning the scattering rate and the energy of a scattered particle. */
    class scatter_process_descriptor
    {
      public:
        scatter_process_descriptor() : initial_energy_(0), final_energy_(0), rate_(0), generation_rate_(0) {}
        //scatter_process_descriptor(double e, double r, double g = 0) : initial_energy_(e), final_energy_(e), rate_(r), generation_rate_(g) {}

        void   initial_energy(double e) { initial_energy_ = e; }
        double initial_energy() const { return initial_energy_; }

        void   final_energy(double e) { final_energy_ = e; }
        double final_energy() const { return final_energy_; }

        void   rate(double r) { rate_ = r; }
        double rate() const { return rate_; }

        void   generation_rate(double r) { generation_rate_ = r; }
        double generation_rate() const { return generation_rate_; }


      private:
        double initial_energy_;
        double final_energy_;
        double rate_;
        double generation_rate_;
    };

    template <typename DeviceType>
    class scattering_base
    {
      protected:
        typedef typename DeviceType::mesh_type              MeshType;

        typedef typename viennagrid::result_of::point<MeshType>::type      PointType;
        typedef typename viennagrid::result_of::vertex<MeshType>::type     VertexType;
        typedef typename viennagrid::result_of::facet<MeshType>::type      FacetType;
        typedef typename viennagrid::result_of::cell<MeshType>::type       CellType;

      public:
        typedef std::vector<scatter_process_descriptor>   scatter_processes_type;
        typedef scatter_processes_type                    value_type;

        scattering_base(DeviceType const & device, viennashe::config const & conf) : device_(device), conf_(conf) {}

        virtual ~scattering_base() {}

        virtual scatter_processes_type operator()(FacetType const & elem,
                                                  double kinetic_energy,
                                                  viennashe::carrier_type_id ctype) const = 0;

        virtual scatter_processes_type operator()(CellType const & elem,
                                                  double kinetic_energy,
                                                  viennashe::carrier_type_id ctype) const = 0;

        virtual scatter_process_id id() const = 0;

      protected:
        DeviceType        const & device_;
        viennashe::config const & conf_;
    };


  } //namespace she
} //namespace viennashe

#endif

