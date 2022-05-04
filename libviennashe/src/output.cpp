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

// C++ includes
#include "libviennashe/src/viennashe_all.hpp"
#include "libviennashe/src/quantity_wrappers.hpp"


#include <iostream>
#include <iosfwd>

// C includes
#include "libviennashe/include/output.h"

namespace libviennashe
{

  /**
   * @brief Writes the given quantity to gnuplot
   * @param sim The ViennaSHE simulator
   * @param reg The quantity register
   * @param name The name of the quantity in the register
   * @param filename The filename (or path) to the gnuplot file (will be overwritten!)
   */
  template < typename SimulatorT >
  void write_to_gnuplot(SimulatorT const & sim, libviennashe::quan_register_internal & reg,  std::string name, std::string filename)
  {
    typedef typename SimulatorT::device_type  DeviceType;
    typedef typename DeviceType::mesh_type    MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type     CellIterator;

    DeviceType const & device = sim.device();

    std::ofstream writer(filename.c_str());

    if (!writer)
    {
      viennashe::log::error() << "write_to_gnuplot(): Cannot open the file '" << filename << "'" << std::endl;
      throw viennashe::io::cannot_open_file_exception(filename);
    }

    writer << "## ViennaSHE - gnuplot output of '" << name << "' in SI units unless otherwise noted" << std::endl;
    writer << "## x_i are the vertex coordinates (meter) and a_i are the quantity values " << std::endl;

    bool first = true;

    const std::size_t dim = static_cast<std::size_t>(PointType::dim);

    CellContainer cells(device.mesh());

    for (CellIterator cit = cells.begin();
         cit != cells.end();
         ++cit)
    {

      std::vector<double> values(1);
      reg.cell_based.get(name).fill_single(std::size_t(cit->id().get()), values);

      if (values.size() == 0) continue;

      // Write preamble
      if (first)
      {
        writer << "# ";
        for (std::size_t i = 0; i < dim; ++i) writer << " x_" << i << " ";
        for (std::size_t i = 0; i < values.size(); ++i)  writer << " a_" << i << " ";
        first = false;
        writer << std::endl;
      }

      // Write values at point
      for (std::size_t i = 0; i < dim; ++i) writer << viennagrid::centroid(*cit)[i] << " ";
      for (std::size_t i = 0; i < values.size(); ++i)  writer << values[i]       << " ";
      writer << std::endl;
    } // for vertices

    viennashe::log::info() << "* write_to_gnuplot(): Writing data to '"
                           << filename
                           << "' (can be viewed with e.g. gnuplot)" << std::endl;

  }

  /**
   * @brief Writes SHE simulator results to a gnuplot file
   * @param sim The SHE simulator
   * @param ctype The carrier type for which to write the results
   * @param filename The name of the gnuplot file (will be overwritten)
   */
  template < typename SimulatorT >
  void write_she_to_gnuplot(SimulatorT const & sim, viennashe::carrier_type_id ctype, std::string filename)
  {
    typedef typename SimulatorT::device_type DeviceType;
    typedef typename DeviceType::mesh_type   MeshType;

    typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
    typedef typename viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

    typedef typename SimulatorT::edf_type                              EDFType;
    typedef typename SimulatorT::generalized_edf_type                  GeneralizedEDFType;
    typedef typename viennashe::config::dispersion_relation_type       DispersionRelationType;

    // Safety check
    if (ctype == viennashe::ELECTRON_TYPE_ID && (!sim.config().with_electrons() || sim.config().get_electron_equation() != viennashe::EQUATION_SHE))
    {
      viennashe::log::warning() << "Warning! write_she_to_gnuplot(): Electron SHE Quantities cannot be written to file. There are no quantities to write!" << std::endl;
      return;
    }
    if (ctype == viennashe::HOLE_TYPE_ID && (!sim.config().with_holes() || sim.config().get_hole_equation() != viennashe::EQUATION_SHE))
    {
      viennashe::log::warning() << "Warning! write_she_to_gnuplot(): Hole SHE Quantities cannot be written to file. There are no quantities to write!" << std::endl;
      return;
    }

    DeviceType        const & device = sim.device();
    viennashe::config const & conf   = sim.config();
    DispersionRelationType    disp   = conf.dispersion_relation(ctype);
    EDFType            edf  = sim.edf(ctype);
    GeneralizedEDFType gedf = sim.generalized_edf(ctype);


    std::ofstream writer(filename.c_str());

    if (!writer)
    {
      throw viennashe::io::cannot_open_file_exception(filename);
    }

    writer << "## ViennaSHE - gnuplot output of SHE "
           << ((ctype == viennashe::ELECTRON_TYPE_ID) ? "Electron" : "Hole")
           <<  " Quantities in SI units unless otherwise noted" << std::endl;
    writer << "## x_i are the vertex coordinates (meter)  " << std::endl;

    bool first = true;

    CellContainer cells(device.mesh());
    for (CellIterator cit = cells.begin();
         cit != cells.end();
         ++cit)
    {
      // Write preamble
      if (first)
      {
        writer << "# ";
        for (std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i) writer << " x_" << i << " ";
        writer << "     energy     edf      generalized_edf      dos      vg" << std::endl;
        first = false;
      }

      // Write values at point
      for (std::size_t index_H = 1; index_H < sim.quantities().carrier_distribution_function(ctype).get_value_H_size()-1; ++index_H)
      {
        const double eps  = sim.quantities().carrier_distribution_function(ctype).get_kinetic_energy(*cit, index_H);
        const double dos  = disp.density_of_states(eps);
        const double velo = disp.velocity(eps);

        if (eps >= 0.0)
        {
          for (std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i) writer << viennagrid::centroid(*cit)[i] << " ";
          writer << eps << " " << edf(*cit, eps, index_H) << " " << gedf(*cit, eps, index_H) << " " << dos << " " << velo << std::endl;
        }
      }

      writer << std::endl;
      writer << std::endl;
    } // for vertices

    viennashe::log::info() << "* write_she_to_gnuplot(): Writing data to '"
      << filename
      << "' (can be viewed with e.g. gnuplot)" << std::endl;
  }

} // namespace libviennashe


#ifdef	__cplusplus
extern "C"
{
#endif

viennasheErrorCode viennashe_write_to_gnuplot(viennashe_quan_register reg, char const * name, char const * filename)
{
  try
  {
    //
    // CHECKS
    if (reg == NULL)
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot: The quantity regiser (reg) must not be NULL!" << std::endl;
      return 1;
    }
    if (filename == NULL)
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot: The filename must not be NULL!" << std::endl;
      return 3;
    }
    // Get the internal structure
    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    // Get internal simulator
    viennashe_simulator_impl const * int_sim = int_reg->int_sim;
    // More checks
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }
    // Check if the quan exists
    if(! int_reg->cell_based.has_quan(name))
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot(): The quantity '" << name << "' does not exist!" << std::endl;
      return 2;
    }

    // Do the actual work
    if(int_sim->stype == libviennashe::meshtype::line_1d)
    {
      libviennashe::write_to_gnuplot(*(int_sim->sim1d), *int_reg, std::string(name), std::string(filename));
    }
    else if(int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      libviennashe::write_to_gnuplot(*(int_sim->simq2d), *int_reg, std::string(name), std::string(filename));
    }
    else if(int_sim->stype == libviennashe::meshtype::triangular_2d)
    {
      libviennashe::write_to_gnuplot(*(int_sim->simt2d), *int_reg, std::string(name), std::string(filename));
    }
    else if(int_sim->stype == libviennashe::meshtype::hexahedral_3d)
    {
      libviennashe::write_to_gnuplot(*(int_sim->simh3d), *int_reg, std::string(name), std::string(filename));
    }
    else if(int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      libviennashe::write_to_gnuplot(*(int_sim->simt3d), *int_reg, std::string(name), std::string(filename));
    }
    else
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot(): Unkown grid type!" << std::endl;
      return -2;
    }
  }
  catch (std::exception const & ex)
  {
    viennashe::log::error() << "ERROR! write_to_gnuplot: Exception!" << std::endl;
    viennashe::log::error() << "What? " << ex.what() << std::endl;
    return -1;
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR! write_to_gnuplot: UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT  viennasheErrorCode viennashe_write_she_results_to_gnuplot(viennashe_quan_register reg, viennashe_carrier_ids carriertype, char const * filename)
{
  try
  {
    //
    // CHECKS
    if (reg == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_write_she_results_to_gnuplot: The quantity regiser (reg) must not be NULL!" << std::endl;
      return 1;
    }
    if (filename == NULL)
    {
      viennashe::log::error() << "ERROR! viennashe_write_she_results_to_gnuplot: The filename must not be NULL!" << std::endl;
      return 3;
    }

    // Get the internal structure
    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    // Get internal simulator
    viennashe_simulator_impl const * int_sim = int_reg->int_sim;
    // More checks
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! viennashe_write_she_results_to_gnuplot(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    viennashe::carrier_type_id ctype = ((carriertype == viennashe_electron_id) ? viennashe::ELECTRON_TYPE_ID : viennashe::HOLE_TYPE_ID);

    // Do the actual work
    if (int_sim->stype == libviennashe::meshtype::line_1d)
    {
      libviennashe::write_she_to_gnuplot(*(int_sim->sim1d), ctype, std::string(filename));
    }
    else if (int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
    {
      libviennashe::write_she_to_gnuplot(*(int_sim->simq2d), ctype, std::string(filename));
    }
    else if (int_sim->stype == libviennashe::meshtype::triangular_2d)
    {
      libviennashe::write_she_to_gnuplot(*(int_sim->simt2d), ctype, std::string(filename));
    }
    else if (int_sim->stype == libviennashe::meshtype::hexahedral_3d)
    {
      libviennashe::write_she_to_gnuplot(*(int_sim->simh3d), ctype, std::string(filename));
    }
    else if (int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
    {
      libviennashe::write_she_to_gnuplot(*(int_sim->simt3d), ctype, std::string(filename));
    }
    else
    {
      viennashe::log::error() << "ERROR! viennashe_write_she_results_to_gnuplot(): Unkown grid type!" << std::endl;
      return -2;
    }
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR! viennashe_write_she_results_to_gnuplot: UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


#ifdef	__cplusplus
}
#endif

