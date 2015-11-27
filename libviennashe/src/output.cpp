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

    viennagrid_dimension geo_dim;
    viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim);
    const std::size_t dim = std::size_t(geo_dim);

    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {

      std::vector<double> values(1);
      reg.cell_based.get(name).fill_single(std::size_t(viennagrid_index_from_element_id(*cit)), values);

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
      viennagrid_numeric centroid[3];
      viennagrid_element_centroid(device.mesh(), *cit, centroid);
      for (std::size_t i = 0; i < dim; ++i) writer << centroid[i] << " ";
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

    viennagrid_dimension geo_dim;
    viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim);
    const std::size_t dim = std::size_t(geo_dim);

    viennagrid_dimension cell_dim;
    viennagrid_mesh_cell_dimension_get(device.mesh(), &cell_dim);

    viennagrid_element_id *cells_begin, *cells_end;
    viennagrid_mesh_elements_get(device.mesh(), cell_dim, &cells_begin, &cells_end);
    for (viennagrid_element_id *cit  = cells_begin;
                                cit != cells_end;
                              ++cit)
    {
      // Write preamble
      if (first)
      {
        writer << "# ";
        for (std::size_t i = 0; i < dim; ++i) writer << " x_" << i << " ";
        writer << "     energy     edf      generalized_edf      dos      vg" << std::endl;
        first = false;
      }

      viennagrid_numeric centroid[3];
      viennagrid_element_centroid(device.mesh(), *cit, centroid);

      // Write values at point
      for (std::size_t index_H = 1; index_H < sim.quantities().carrier_distribution_function(ctype).get_value_H_size()-1; ++index_H)
      {
        const double eps  = sim.quantities().carrier_distribution_function(ctype).get_kinetic_energy(*cit, index_H);
        const double dos  = disp.density_of_states(eps);
        const double velo = disp.velocity(eps);

        if (eps >= 0.0)
        {
          for (std::size_t i = 0; i < dim; ++i) writer << centroid[i] << " ";
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

    // Check if the quan exists
    if(! int_reg->cell_based.has_quan(name))
    {
      viennashe::log::error() << "ERROR! write_to_gnuplot(): The quantity '" << name << "' does not exist!" << std::endl;
      return 2;
    }

    libviennashe::write_to_gnuplot(int_sim->sim_, *int_reg, std::string(name), std::string(filename));
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

    viennashe::carrier_type_id ctype = ((carriertype == viennashe_electron_id) ? viennashe::ELECTRON_TYPE_ID : viennashe::HOLE_TYPE_ID);

    libviennashe::write_she_to_gnuplot(int_sim->sim_, ctype, std::string(filename));
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

