#ifndef VIENNASHE_MATERIALS_MATERIALS_HPP
#define VIENNASHE_MATERIALS_MATERIALS_HPP

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
#include <stdexcept>
#include <string>
#include <algorithm>

// viennashe::materials
#include "viennashe/materials/exception.hpp"

// viennashe
#include "viennashe/forwards.h"

// viennashe::physics
#include "viennashe/physics/constants.hpp"


/** @file viennashe/materials/all.hpp
    @brief A very simple material database. Needs to be replaced by something more versatile soon
*/

namespace viennashe
{
  /** @brief Namespace containing a small materials library. */
  namespace materials
  {
    /** @brief A class referring to any metal contact. */
    struct metal
    {
      enum { id = 100 };
    };

    /** @brief A class referring to silicon and providing certain material parameters Note that this is the default material! */
    struct si
    {
      enum { id = 0 };  //Note that this is the default material!

      static double permittivity() { return 11.9 * viennashe::physics::constants::eps_0; }

      //
      // carrier type agnostic versions:
      //
      static double dos_effective_mass(viennashe::carrier_type_id ctype)
      {
        if (ctype == ELECTRON_TYPE_ID)
          return 0.32 * viennashe::physics::constants::mass_electron;
        else
          return 0.8 * viennashe::physics::constants::mass_electron;
      }

      static double longitudinal_effective_mass(viennashe::carrier_type_id ctype)
      {
        if (ctype == ELECTRON_TYPE_ID)
          return 0.98 * viennashe::physics::constants::mass_electron;
        else
          return 0.16 * viennashe::physics::constants::mass_electron;
      }

      static double transverse_effective_mass(viennashe::carrier_type_id ctype)
      {
        if (ctype == ELECTRON_TYPE_ID)
          return 0.19 * viennashe::physics::constants::mass_electron;
        else
          return 0.16 * viennashe::physics::constants::mass_electron;
      }

      static double conductivity_effective_mass(viennashe::carrier_type_id ctype)
      {
        if (ctype == ELECTRON_TYPE_ID)
        {
          return 3.0 / (  1.0/longitudinal_effective_mass(ctype)
                        + 1.0/transverse_effective_mass(ctype)
                        + 1.0/transverse_effective_mass(ctype) );
        }
        else
        {
          //      return 3.0 / (  1.0/longitudinal_effective_mass(viennashe::hole_tag())
          //                    + 1.0/transverse_effective_mass(viennashe::hole_tag())
          //                    + 1.0/transverse_effective_mass(viennashe::hole_tag()) );
          return 0.36 * viennashe::physics::constants::mass_electron;         //Note: This should be the harmonic mean for consistency...
        }
      }

      static double band_gap() { return 1.11 * viennashe::physics::constants::q; /* in Joule @ 300 K */ }

      static double get_mass_density() { return 2330.0; /* kg/m^3 */ }

      static double specific_heat() { return 1.6e6; /* J m^-3 K^-1 */ }
      static double reference_temperature() { return 300.0; /* K */ }

      static double diffusivity() { return 10; /* TODO: real values */ }
    };

    /** @brief A class referring to silicon dioxide */
    struct sio2
    {
      enum { id = 1 };
      static double permittivity() { return 3.9 * viennashe::physics::constants::eps_0; }
      static double diffusivity() { return 100; /* TODO: real values */ }
    };

    /** @brief A class referring to hafnium dioxide */
    struct hfo2
    {
      enum { id = 2 };
      static double permittivity() { return 25.0 * viennashe::physics::constants::eps_0; }
      static double diffusivity() { return 100; /* TODO: real values */ }
    };

    /** @brief A class referring to silicon oxynitride  */
    struct sion
    {
      // Find a way to incorporate nitrogen and oxygen content. The relative permitivity of SiO_xN_y depends on x and y!
      enum { id = 3 };
      static double permittivity() { return 10.0 * viennashe::physics::constants::eps_0; }
      static double diffusivity() { return 100; /* TODO: real values */ }
    };

    /** @brief A class referring to an ideal gas */
    struct gas
    {
      enum { id = 4 };
      static double permittivity() { return 1.0 * viennashe::physics::constants::eps_0; }
      static double diffusivity() { return 100; /* TODO: real values */ }
    };

    //
    // Convenience functions:
    //

    /** @brief Convenience function for checking whether the supplied material ID refers to a metal */
    inline bool is_conductor(long material_id)
    {
      if (material_id == metal::id)
        return true;

      return false;
    }

    /** @brief Convenience function for checking whether the supplied material ID refers to a semiconductor */
    inline bool is_semiconductor(long material_id)
    {
      if (material_id == si::id)
        return true;

      return false;
    }

    /** @brief Convenience function for checking whether the supplied material ID refers to an oxide */
    inline bool is_insulator(long material_id)
    {
      if (material_id == sio2::id)
        return true;
      if (material_id == hfo2::id)
        return true;
      if (material_id == sion::id)
        return true;
      if (material_id == gas::id)
        return true;

      return false;
    }

    /** @brief Simple checker class for checking whether a certain material is of a given type (conductor, semiconductor, insulator).
      *
      * Uses a runtime dispatch, since this checker is not considered to be performance critical, yet provides the benefit of smaller compilation times.
      */
    class checker
    {
    public:
      checker(material_category_id category) : category_(category) {}

      bool operator()(long material_id) const
      {
        switch (category_)
        {
        case MATERIAL_CONDUCTOR_ID:        return  is_conductor(material_id);
        case MATERIAL_NO_CONDUCTOR_ID:     return !is_conductor(material_id);
        case MATERIAL_SEMICONDUCTOR_ID:    return  is_semiconductor(material_id);
        case MATERIAL_NO_SEMICONDUCTOR_ID: return !is_semiconductor(material_id);
        case MATERIAL_INSULATOR_ID:        return  is_insulator(material_id);
        default:                           return !is_insulator(material_id); //only remaining option
        }

        //std::stringstream ss;
        //ss << "Invalid material category encountered!" << std::endl;
        //throw invalid_material_exception(ss.str());
      }
    private:
      material_category_id category_;
    };


    //
    // Permittivity accessor:
    //

    /** @brief Convenience function for returning the permittivity of the material identified by the ID provided */
    inline double permittivity(long material_id)
    {
      switch (material_id)
      {
        case si::id:
           return si::permittivity();
        case sio2::id:
           return sio2::permittivity();
        case hfo2::id:
           return hfo2::permittivity();
        case sion::id:
           return sion::permittivity();
        case gas::id:
           return gas::permittivity();
      };

      std::stringstream ss;
      ss << "Cannot retrieve permittivity from material with id: " << material_id << std::endl;
      throw invalid_material_exception(ss.str());
    }

    inline double diffusivity(long material_id)
    {
      switch (material_id)
      {
        case si::id:
           return si::diffusivity();
        case sio2::id:
           return sio2::diffusivity();
        case hfo2::id:
           return hfo2::diffusivity();
        case sion::id:
           return sion::diffusivity();
      };

      std::stringstream ss;
      ss << "Cannot retrieve diffusivity from material with id: " << material_id << std::endl;
      throw invalid_material_exception(ss.str());
    }


    /** @brief Returns the material ID for a material identified by a string
     *
     * @param material_name The name of the material. Needs to be by value, because of std::transform
     */
    inline long get_material_id(std::string material_name)
    {
      // to upper to ease the comparison and to be case insensitive
      std::transform(material_name.begin(), material_name.end(), material_name.begin(), ::toupper);

      if (material_name == "SI")
      {
          return si::id;
      }
      else if (material_name == "SIO2")
      {
          return sio2::id;
      }
      else if (material_name == "HFO2")
      {
          return hfo2::id;
      }
      else if (material_name == "SION")
      {
          return sion::id;
      }
      else if (material_name == "GAS")
      {
          return gas::id;
      }
      else if (material_name == "METAL")
      {
          return metal::id;
      }
      else
      {
        std::stringstream ss;
        ss << "Cannot retrieve material ID for material named: '" << material_name << "'" << std::endl;
        throw invalid_material_exception(ss.str());
      }
    }

  } //namespace materials
} //namespace viennashe

#endif

