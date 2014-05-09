#ifndef VIENNASHE_TESTS_MOS1D_DG_HPP
#define	VIENNASHE_TESTS_MOS1D_DG_HPP
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

/** @file tests/src/mos1d_dg.hpp
    @brief Contains common code for the 1D MOS DG tests
 */

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include <iostream>
#include <cstdlib>
#include <vector>

#include "tests/src/common.hpp"

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default configuration and centroid() algorithm:
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/algorithm/centroid.hpp"


/** \file mos1d_dg.hpp This is a header file containing all common definitions and
 *                     declarations for the density gradient model tests on 1D MOS structures.
 */


/**
 * @brief MOS 1D test-grid generator
 * @param len_gate The length of the Gate in meter
 * @param cs_gate Spacing of two cell-centroids in meter
 * @param len_oxide The length of the Oxide in meter
 * @param cs_ox Grid spacing in the Oxide
 * @param len_bulk Lenght of the Bulk in meter
 * @param cs_bulk Grid spacing for the Bulk
 */
struct mos1d_mesh_generator
{

  mos1d_mesh_generator(double len_gate, double cs_gate, double len_oxide, double cs_ox, double len_bulk, double cs_bulk)
    : len_gate_(len_gate), cs_gate_(cs_gate), len_oxide_(len_oxide), cs_ox_(cs_ox), len_bulk_(len_bulk), cs_bulk_(cs_bulk)
  { }

  template < typename MeshT, typename SegmentationT >
  void operator()(MeshT & mesh, SegmentationT & seg) const
  {
    viennashe::util::device_generation_config gconf;

    gconf.add_segment(0,                              len_gate_,   static_cast<unsigned long>(std::ceil(len_gate_ / cs_gate_)) + 1);
    gconf.add_segment(len_gate_,                      len_oxide_,  static_cast<unsigned long>(std::ceil(len_oxide_ / cs_ox_) )    );
    gconf.add_segment(len_gate_+len_oxide_,           +20e-9,      static_cast<unsigned long>(std::ceil(+20e-9 / 0.05e-9)    )    );
    gconf.add_segment(len_gate_+len_oxide_+20e-9,     len_bulk_,   static_cast<unsigned long>(std::ceil(len_bulk_ / cs_bulk_))    );
    gconf.add_segment(len_gate_+len_oxide_+20e-9+len_bulk_,  2e-9, static_cast<unsigned long>(std::ceil(2e-9 / cs_bulk_)     ) + 1);

    viennashe::util::generate_device(mesh, seg, gconf);
  }

private:
  double len_gate_;
  double cs_gate_;
  double len_oxide_;
  double cs_ox_;
  double len_bulk_;
  double cs_bulk_;

};


/** @brief Initalizes the MOS 1D device for testing.
*
* @param device The device class that is to be initalized
* @param Vg_init The initial gate voltage
* @param Nd The donor doping (m^-3)
* @param Na The acceptor doping (m^-3)
*/
template <typename DeviceType>
void init_device(DeviceType & device, double Vg_init, double Nd, double Na)
{
  typedef typename DeviceType::segment_type        SegmentType;

  SegmentType const & gate     = device.segment(0);
  SegmentType const & oxide    = device.segment(1);
  SegmentType const & silicon  = device.segment(2);
  SegmentType const & silicon2 = device.segment(3);
  SegmentType const & bulk     = device.segment(4);

  std::cout << "* init_device(): Setting material ..." << std::endl;

  device.set_material(viennashe::materials::metal(), gate);
  device.set_material(viennashe::materials::metal(), bulk);
  device.set_material(viennashe::materials::sio2(),  oxide);
  device.set_material(viennashe::materials::si(),    silicon);
  device.set_material(viennashe::materials::si(),    silicon2);

  std::cout << "* init_device(): Setting doping (per cell) ..." << std::endl;

  device.set_doping_n(Nd, silicon);
  device.set_doping_p(Na, silicon);
  device.set_doping_n(Nd, silicon2);
  device.set_doping_p(Na, silicon2);

  std::cout << "* init_device(): Setting contact potentials (per cell) ..." << std::endl;

  device.set_contact_potential(0.0,      bulk);
  device.set_contact_potential(Vg_init,  gate);

  std::cout << "* init_device(): DONE!" << std::endl;
} // init_device()

/** @brief A simple reference data reader */
struct reference_values
{
  typedef std::map<long, std::vector<double> > storage_type;

  /**
   * @brief Reads reference data from a CSV text file, where the main index to the data is the cell-id
   * @param filename The path/filename of the input file (text file)
   * @param delimiter The value delimiter in the text file, usually it's just a space character
   * @return True on success else false
   */
  bool read_file(std::string filename, char delimiter = ' ')
  {
    std::ifstream       reader(filename.c_str());
    if(!reader) return false;

    long i = 0;
    while(reader.good() && !reader.eof())
    {
      std::vector<double> result;
      std::string         line;
      std::getline(reader,line);

      std::stringstream          lineStream(line);
      std::string                cell;

      if(line.length() == 0) continue;
      if(line[0] == '#' ) continue;

      while(std::getline(lineStream, cell, delimiter))
      {
        result.push_back(atof(cell.c_str()));
      }

      data_.insert(std::make_pair(i,result));
      i++;
    }
    return true;
  }

  // Returns a map to the reference data ... index is the cell-id
  storage_type const & data() const { return data_; }

private:
  storage_type data_;
};

/**
 * @brief Tests the given quantity against reference data in a file. Uses reference_values
 * @param quan An accessor (cell-based) to the quantity
 * @param device Reference to the device
 * @param filename Path/Filename of the file (CSV, sperator is space) with reference values
 * @return True if the given data (quan) matches the reference from file (filename)
 */
template <typename DeviceType, typename AccessorType >
bool test_result(AccessorType const & quan, DeviceType const & device, std::string filename )
{
  typedef typename DeviceType::mesh_type           MeshType;

  typedef typename viennagrid::result_of::point<MeshType>::type     PointType;

  typedef typename viennagrid::result_of::const_cell_range<MeshType>::type   CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

  const double tol = 1e-1;

  std::cout << "test_result(): file = '" << filename << "'" << std::endl;
  try
  {
    reference_values referencefile;

    bool ok = referencefile.read_file(filename);
    if(!ok)
    {
      std::cerr << "test_result(): Unable to read file '" << filename << "'" << std::endl;
      return false;
    }

    long linenum = 0;
    CellContainer cells(device.mesh());
    for (CellIterator cit = cells.begin();
         cit != cells.end();
         ++cit, ++linenum )
    {
      //write values at point
      for (std::size_t i = 0; i < static_cast<std::size_t>(PointType::dim); ++i)
      {
        const double value    = viennagrid::centroid(*cit)[i];
        const double refvalue = referencefile.data().at(linenum)[i];
        if(!viennashe::testing::fuzzy_equal(refvalue, value, tol))
        {
          std::cerr << "test_result(): Failure at " << *cit << " => " << quan(*cit) << std::endl;
          return false;
        }
      }
      const double value    = quan(*cit);
      const double refvalue = referencefile.data().at(linenum)[PointType::dim];
      if(!viennashe::testing::fuzzy_equal(refvalue, value, tol))
      {
        std::cerr << "test_result(): Failure at " << *cit << " => " << quan(*cit) << std::endl;
        return false;
      }
    } // for vertices
  }
  catch(std::exception & ex)
  {
    std::cerr << "test_result(): FAILED! What? " << ex.what() << std::endl;
    return false;
  }

  return true;
}



#endif	/* VIENNASHE_TESTS_MOS1D_DG_HPP */

