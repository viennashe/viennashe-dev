#ifndef VIENNASHE_MODELS_HOT_CARRIER_HPP
#define VIENNASHE_MODELS_HOT_CARRIER_HPP

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
#include "viennashe/models/exception.hpp"
#include "viennashe/physics/constants.hpp"
#include "viennashe/config.hpp"
#include "viennashe/math/integrator.hpp"

#include "viennashe/postproc/electric_field.hpp"
#include "viennashe/simulator.hpp"

/** @file viennashe/models/hot_carrier.hpp
    @brief Implements hot carrier model by Tyaginov et al. (cf. Tyaginov et al. IRW 2013)
*/

namespace viennashe
{
  namespace models
  {
    /** @brief Namespace to hide implementation details */
    namespace detail
    {
      /**
       * @brief Encapsulates the definition of the Keldysh capture cross section: ( (E - Eadd) / Eref )^p * sigma0
       * @param Eadd Energy in Joule: ( (E - Eadd) / Eref )^p * sigma0
       * @param Eref Energy in Joule: ( (E - Eadd) / Eref )^p * sigma0
       * @param sigma_0 Constant capture cross section: ( (E - Eadd) / Eref )^p * sigma0
       * @param p Dimensionless exponent
       */
      struct capture_cross_section_hcd
      {
      public:
        capture_cross_section_hcd(double Eadd, double Eref, double sigma_0, double p = 1)
            : Eadd_(Eadd), Eref_(Eref), sigma0_(sigma_0), p_(p)
        {
          if(Eref_ <= 0.0) Eref_ = 1.0;
          if(p_ <= 0.0) p_ = 1.0;
        }

        double operator()(double E) const { return sigma0_ * std::pow( (E + Eadd_)/Eref_ , p_);  }

        void set_Eadd(double E) { this->Eadd_ = E; }

      private:
        double Eadd_;
        double Eref_;
        double sigma0_;
        double p_;

      };

      template < typename EDFAccessorT, typename ElementType >
      class acceleration_integrand
      {
        public:
          typedef viennashe::config::dispersion_relation_type      dispersion_relation_type;

          acceleration_integrand(dispersion_relation_type const & dispersion,
                                 EDFAccessorT const & edf,
                                 ElementType const & el,
                                 capture_cross_section_hcd const & ccs)
            : dispersion_(dispersion), edf_(edf), el_(el), ccs_(ccs) {}

          double operator()(double kinetic_energy) const
          {
            /*
            const double ccs = ccs_(kinetic_energy);
            const double dos = dispersion_.density_of_states(kinetic_energy);
            const double vg  = dispersion_.velocity(kinetic_energy);
            const double edf = edf_(el_, kinetic_energy);
             */
            return (  std::max(ccs_(kinetic_energy), 0.0)
                    * dispersion_.symmetry_factor()     //in order to get consistent results when dividing by the carrier density
                    * dispersion_.density_of_states(kinetic_energy)
                    * dispersion_.velocity(kinetic_energy)
                    * std::abs(edf_(el_, kinetic_energy))
                   );
          }

        private:
          dispersion_relation_type const & dispersion_;
          EDFAccessorT const & edf_;
          ElementType const & el_;
          capture_cross_section_hcd const & ccs_;
      };
    } // namespace detail

    /**
     * @brief Plain-Old-Data type holding all parameters for the hot carrier model.
     */
    struct hcd_parameters
    {
      double Ea_mean;
      double Ea_sigma;
      double sigma0;
      double Ep;
      double p;
      double Eref;
      double omega;
      double hbaromega;
      double nu;
      double Nitmax;
      double N;
      double Ethn, Ethp;

      hcd_parameters()
        : Ea_mean(viennashe::physics::constants::q * 1.5),
          Ea_sigma(viennashe::physics::constants::q * 0.25),
          sigma0(1e-22),
          Ep(viennashe::physics::constants::q * 0.2),
          p(11),
          Eref(viennashe::physics::constants::q),
          omega(1e11),
          hbaromega(viennashe::physics::constants::q * 0.075),
          nu(1e-10),//(1e-13),
          Nitmax(1e19),
          N(3),
          Ethn(viennashe::physics::constants::q * 0.1),
          Ethp(viennashe::physics::constants::q * 0.1)
      {}
    };

    template < typename DeviceT, typename SHEQuanT >
    struct acceleration_integral
    {
    private:
      typedef typename DeviceT::mesh_type           MeshType;

      typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type   VertexContainer;
      typedef typename viennagrid::result_of::const_vertex_range<MeshType>::type   CellContainer;

    public:

      typedef typename viennagrid::result_of::vertex<MeshType>::type       VertexType;
      typedef typename viennagrid::result_of::edge<MeshType>::type         EdgeType;
      typedef typename viennagrid::result_of::cell<MeshType>::type         CellType;

      typedef typename viennashe::she::edf_wrapper<DeviceT, SHEQuanT> edf_type;

      acceleration_integral(viennashe::config const & conf,
                            SHEQuanT const & quan)
        : conf_(conf), quan_(quan), edf_(conf, quan) { }

      template < typename CaptureCrossSectionT >
      double operator()( CellType const & elem, CaptureCrossSectionT const & ccs, double Eth) const
      {
        double     res = 0;
        typename viennashe::config::dispersion_relation_type dispersion = conf_.dispersion_relation(quan_.get_carrier_type_id());
        detail::acceleration_integrand<edf_type, CellType> acc_int(dispersion, edf_, elem, ccs);
        for (std::size_t index_H=1; index_H < quan_.get_value_H_size() - 1; ++index_H)
        {
          const double energy_mid = quan_.get_kinetic_energy(elem, index_H);
          const double energy_lower = std::max( (quan_.get_kinetic_energy(elem, index_H - 1) + energy_mid) / 2.0, 0.0);
          const double energy_upper =           (quan_.get_kinetic_energy(elem, index_H + 1) + energy_mid) / 2.0;

          if (energy_upper < 0 || energy_lower < Eth)
            continue;

          if (energy_lower < energy_upper)
          {
            res += viennashe::math::integrate(acc_int, energy_lower, energy_upper, viennashe::math::IntGauss<5>());
          }
        }
        return  res;
      }

    private:

      viennashe::config conf_;
      SHEQuanT quan_;
      edf_type edf_;
    };



    template < typename DeviceType, typename AccelerationIntegralElectronsT, typename AccelerationIntegralHolesT >
    class hcd_model
    {
    private:

      typedef typename DeviceType::mesh_type           MeshType;

    public:
      typedef typename viennagrid::result_of::cell<MeshType>::type         CellType;

      typedef typename viennagrid::result_of::point<MeshType>::type    PointType;

      typedef hcd_model<DeviceType, AccelerationIntegralElectronsT, AccelerationIntegralHolesT>    self_type;
      typedef DeviceType device_type;


      hcd_model(DeviceType const & device,
                AccelerationIntegralElectronsT const & for_electrons,
                AccelerationIntegralHolesT const & for_holes,
                hcd_parameters const & params)
       : device_(device), acc_electrons_(for_electrons), acc_holes_(for_holes),  params_(params) { }

      hcd_model(hcd_model const & o) : device_(o.device_), acc_electrons_(o.acc_electrons_),
                                       acc_holes_(o.acc_holes_), params_(o.params_)
      {  }

      template < typename ElectricFieldAccessorT >
      double operator()(CellType const & cell, double time, PointType dipole, ElectricFieldAccessorT const & Eox) const
      {
        return this->evaluate(cell, time, dipole, Eox);
      }


      template < typename ElectricFieldAccessorT >
      double evaluate(CellType const & cell, double time,
                      PointType dipole,
                      ElectricFieldAccessorT const & Eox) const
      {
        typedef std::map<double, double>          energy_map_type;
        typedef energy_map_type::const_iterator   energy_map_iterator;

        if (time < 0) throw viennashe::models::invalid_parameter_exception("Time must be greater or equal to zero!", time);

        double Nt = 0;

        energy_map_type activation_energies;
        this->sample_activation_energy(activation_energies, 100, params_.Ea_mean, params_.Ea_sigma, params_.Nitmax);

        const double dipole_reduction = this->get_dipole_reduction(cell, dipole, Eox);

        detail::capture_cross_section_hcd ccs_mve( - params_.hbaromega, params_.Eref, params_.sigma0, 1);
        detail::capture_cross_section_hcd ccs_ab( 0, params_.Eref, params_.sigma0, params_.p);

        const double kup   = this->get_kup(cell, ccs_mve);
        const double kdown = this->get_kdown(cell, ccs_mve);

        // for each activation energy sample
        for (energy_map_iterator it = activation_energies.begin();
             it != activation_energies.end();
             ++it)
        {
          const double Ea     = it->first;
          const double Nitmax = it->second;

          double C = 0;
          double R = 0;
          double P = 0;

          // for each level of the harmonic oscillator
          for (std::size_t i = 0; i < params_.N; ++i)
          {
            const double r = this->get_r(i, Ea, dipole_reduction, cell, ccs_ab);
            //if (r < 0)  //TODO: For some reason unknown to KR this can be negative. Skipping such cases for now...
            //  continue;

            const double p = this->get_p(i, cell);

            P += p;
            C += 1.0 * std::pow(kdown / kup, static_cast<double>(i));
            R += r   * std::pow(kup / kdown, static_cast<double>(i));
          }
          R = R * C; // normalization

          const double hlp = std::sqrt( R * R + 4.0 * P * R * Nitmax);
          if (hlp != hlp)
          {
            std::cerr << "hlp is NaN!" << std::endl;
            std::cerr << "R: " << R << std::endl;
            std::cerr << "P: " << P << std::endl;
            std::cerr << "Nitmax: " << Nitmax << std::endl;

            throw std::runtime_error("hlp is NaN!");
          }
          //std::cout << "hlp: " << hlp << ", P: " << P << std::endl;

          Nt += std::max(( hlp * std::tanh(0.5 * time * hlp) - R ) / (2.0 * P), 0.0);
        }
        //std::cout << "OUTER LOOP TOOK: " << timer.elapsed() << std::endl;

        // TODO: Nt can get negative !!!
        return (Nt > 0 ? Nt : 0);
      }

    protected:

      void sample_activation_energy(std::map<double, double> & activation_energies,
                                    std::size_t num_samples,
                                    double mean, double sigma, double Nitmax) const
      {
        const double C = 1.0/std::sqrt(2*viennashe::math::constants::pi);
        const double energy_step = (sigma * 5.0 * 2.0)/(num_samples);
        double E = mean - energy_step * num_samples * 0.5;
        double norm = 0;
        // Get energies and unnormalized weight factor
        for (std::size_t i = 0; i < num_samples; ++i)
        {
          const double dE = E - mean;
          const double val = C/sigma * exp(-dE*dE/(2*sigma*sigma));

          activation_energies.insert(std::make_pair(E, val));
          norm += val;

          E += energy_step;
        }
        // normalize and scale weight factors = represented max(Nit)
        for(std::map<double,double>::iterator it = activation_energies.begin();
            it != activation_energies.end(); ++it)
        {
          it->second = it->second / norm * Nitmax;
        }

      }

      double get_r(std::size_t i, double Ea, double dipole_energy, CellType const & cell, detail::capture_cross_section_hcd & ccs_ab) const
      {
        const double kbtl = this->device_.get_lattice_temperature(cell) * viennashe::physics::constants::kB;
        const double Ei   = i * params_.hbaromega;
        const double Eadd = Ea - dipole_energy - Ei;

        //detail::capture_cross_section_hcd ccs_ab( -1.0 * Eadd, params_.Eref, params_.sigma0, params_.p);
        ccs_ab.set_Eadd(-1.0 * Eadd);

        const double i_ab_electrons = acc_electrons_(cell, ccs_ab, params_.Ethn);
        const double i_ab_holes = acc_holes_(cell, ccs_ab, params_.Ethp);
        const double i_ab_value = i_ab_electrons + i_ab_holes;

        if (i_ab_electrons < 0)
          std::cerr << "i_ab_electrons: " << i_ab_electrons << std::endl;

        if (i_ab_holes < 0)
          std::cerr << "i_ab_holes: " << i_ab_holes << std::endl;

        return i_ab_value + params_.nu * std::exp( - Eadd / kbtl );
      }

      double get_p(std::size_t i, CellType const & cell) const
      {
        if ( i == params_.N - 1)
        {
          const double kbtl = this->device_.get_lattice_temperature(cell) * viennashe::physics::constants::kB;
          return params_.nu * std::exp( - params_.Ep / kbtl);
        }
        return 0;
      }

      double get_kdown(CellType const & cell, detail::capture_cross_section_hcd const & ccs_mve) const
      {
        //detail::capture_cross_section_hcd ccs_mve( - params_.hbaromega, params_.Eref, params_.sigma0, 1);

        const double i_mve_value = acc_electrons_(cell, ccs_mve, params_.Ethn) + acc_holes_(cell, ccs_mve, params_.Ethp) ;

        return i_mve_value + params_.omega ;
      }

      double get_kup(CellType const & cell, detail::capture_cross_section_hcd const & ccs_mve) const
      {
        const double kbtl = this->device_.get_lattice_temperature(cell) * viennashe::physics::constants::kB;

        //detail::capture_cross_section_hcd ccs_mve( - params_.hbaromega, params_.Eref, params_.sigma0, 1);

        const double i_mve_value = acc_electrons_(cell, ccs_mve, params_.Ethn) + acc_holes_(cell, ccs_mve, params_.Ethp) ;

        return i_mve_value + params_.omega * std::exp( - params_.hbaromega / kbtl );
      }

      template < typename ElectricFieldAccessorT >
      double get_dipole_reduction(CellType const & cell, PointType const & direction, ElectricFieldAccessorT const & Eox) const
      {
        typedef typename ElectricFieldAccessorT::value_type accessor_vector_type;

        double redux = 0;
        accessor_vector_type field = Eox(cell);

        for (std::size_t i = 0; i < direction.size(); ++i)
        {
          redux += direction[i] * field[i];
        }
        return redux;
      }

    private:
      DeviceType const & device_;
      AccelerationIntegralElectronsT acc_electrons_;
      AccelerationIntegralHolesT     acc_holes_;
      hcd_parameters params_;
    };

    template < typename SimulatorT >
    class hcd_model_operator
    {
    private:
      typedef typename SimulatorT::device_type::mesh_type   MeshType;

    public:

      typedef typename SimulatorT::device_type        device_type;
      typedef typename SimulatorT::she_quantity_type  she_quantity_type;
      typedef typename SimulatorT::potential_type     potential_type;

      typedef typename viennashe::models::acceleration_integral<device_type, she_quantity_type>  acceleration_integral_type;
      typedef typename viennashe::models::hcd_model<device_type, acceleration_integral_type,
                                                    acceleration_integral_type >                 hcd_model_type;

      typedef typename viennagrid::result_of::point<MeshType>::type    PointType;

      hcd_model_operator(SimulatorT const & she_simulator,
                         hcd_parameters params)
        : she_simulator_(she_simulator),
          params_(params),
          acc_n_(she_simulator.config(), she_simulator.quantities().electron_distribution_function()),
          acc_p_(she_simulator.config(), she_simulator.quantities().hole_distribution_function()),
          hcd_(she_simulator.device(), acc_n_, acc_p_, params_),
          Efield_(she_simulator.device(), she_simulator_.potential())
      {}

      template < typename CellT >
      double operator()(CellT const & cell, PointType dipole, double time) const
      {
        return hcd_.evaluate(cell, time, dipole, Efield_);
      }

    private:

      SimulatorT const & she_simulator_;
      hcd_parameters params_;
      acceleration_integral_type acc_n_;
      acceleration_integral_type acc_p_;
      hcd_model_type             hcd_;
      viennashe::electric_field_wrapper<device_type, potential_type> Efield_; //(device, she_simulator_.potential());
    };


  } // namespace models
} // namespace viennashe


#endif /* VIENNASHE_MODELS_HOT_CARRIER_HPP */

