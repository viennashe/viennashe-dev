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

#if defined(_MSC_VER)
  // Disable name truncation warning obtained in Visual Studio
  #pragma warning(disable:4503)
#endif

#include <iostream>
#include <cstdlib>
#include <vector>

// ViennaSHE includes:
#include "viennashe/core.hpp"

// ViennaGrid default configurations:
#include "viennagrid/config/default_configs.hpp"


/** @brief A field containter for values of ValueType. Supports the accessor interface */
template <typename ValueType>
struct field_container
{
  typedef ValueType   value_type;

  field_container(std::size_t size, value_type default_value = value_type()) : data_(size, default_value) {}

  template <typename ElementType>
  value_type operator()(ElementType const & elem) const
  {
    return data_.at(std::size_t(elem.id().get()));
  }

  template <typename ElementType>
  void operator()(ElementType const & elem, value_type val)
  {
    data_.at(std::size_t(elem.id().get())) = val;
  }

  std::vector<value_type> const & container() const { return data_; }

  std::vector<value_type> data_;
};

/** @brief Tests the quantity transfer between elements and writes the result to output_filename
 * @return True if the test succeeded
 */
template <typename DeviceType>
int test(DeviceType const & device, std::string output_filename)
{
  typedef typename DeviceType::mesh_type     MeshType;

  typedef typename viennagrid::result_of::point<MeshType>::type                 PointType;
  typedef typename viennagrid::result_of::facet<MeshType>::type                 FacetType;
  typedef typename viennagrid::result_of::cell<MeshType>::type                  CellType;

  typedef typename viennagrid::result_of::const_facet_range<MeshType>::type     FacetContainer;
  typedef typename viennagrid::result_of::iterator<FacetContainer>::type        FacetIterator;
  typedef typename viennagrid::result_of::const_cell_range<MeshType>::type      CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type         CellIterator;

  typedef typename viennagrid::result_of::const_facet_range<CellType>::type     FacetOnCellContainer;

  std::vector<double> reference_field(3);
  reference_field[0] = 1.0;
  reference_field[1] = 2.0;
  reference_field[2] = 3.0;

  MeshType const & mesh = device.mesh();

  //
  // Write normal components of field (1, 1, 0) to each edge:
  //
  std::cout << "Writing normal projections on box interfaces to facets..." << std::endl;
  FacetContainer facets(mesh);
  field_container<double> field_component_on_edge_container(facets.size());

  for (FacetIterator fit  = facets.begin();
                     fit != facets.end();
                   ++fit)
  {
    typedef typename viennagrid::result_of::const_coboundary_range<MeshType, FacetType, CellType>::type    CellOnFacetContainer;
    //typedef typename viennagrid::result_of::iterator<CellOnFacetContainer>::type                           CellOnFacetIterator;

    CellOnFacetContainer cells_on_facet(mesh, viennagrid::handle(mesh, *fit));

    // Reference orientation of facet is from first coboundary cell to second coboundary cell:

    PointType n = viennashe::util::outer_cell_normal_at_facet(cells_on_facet[0], *fit);
    n /= viennagrid::norm_2(n); // n should be normalized already, but let's ensure this is indeed the case

    double field_in_n = reference_field[0] * n[0];
    if (PointType::dim >= 2)
      field_in_n += reference_field[1] * n[1];
    if (PointType::dim >= 3)
      field_in_n += reference_field[2] * n[2];
    field_component_on_edge_container(*fit, field_in_n);
  }


  //
  // Now transfer from normal projection on edges to vertices:
  //
  std::cout << "Transferring from edges to vertices..." << std::endl;
  CellContainer cells(mesh);
  field_container<std::vector<double> > my_cell_field_container(cells.size(), std::vector<double>(PointType::dim));

  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    FacetOnCellContainer facets_on_cell(*cit);

    viennashe::util::dual_box_flux_to_cell(device,
                                           *cit,                     facets_on_cell,
                                           my_cell_field_container,  field_component_on_edge_container);

    std::vector<double> e_field(3);
    e_field[0] = my_cell_field_container(*cit)[0];
    if (PointType::dim >= 2)
      e_field[1] = my_cell_field_container(*cit)[1];
    if (PointType::dim >= 3)
      e_field[2] = my_cell_field_container(*cit)[2];

    //
    // Run checks:
    //
    double diff  = (e_field[0] - reference_field[0]) * (e_field[0] - reference_field[0]);
    if (PointType::dim >= 2)
           diff += (e_field[1] - reference_field[1]) * (e_field[1] - reference_field[1]);
    if (PointType::dim >= 3)
           diff += (e_field[2] - reference_field[2]) * (e_field[2] - reference_field[2]);
    if (std::abs(diff) > 1e-10)
    {
      std::cerr << "Test failed at cell " << *cit << std::endl;
      std::cerr << "Diff: " << diff << ", computed vector: (" << e_field[0] << ", " << e_field[1] << ", " << e_field[2] << ")" << std::endl;
      return EXIT_FAILURE;
    }
  }


  //
  // Write to VTK
  //

  std::cout << "Writing to VTK..." << std::endl;
  viennagrid::io::vtk_writer<MeshType> my_vtk_writer;
  my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(my_cell_field_container.container()), "flux");
  my_vtk_writer(mesh, output_filename);

  //
  // Test result:
  //

  return EXIT_SUCCESS;
}



int main()
{
  typedef viennagrid::line_1d_mesh           MeshType1d;
  typedef viennashe::device<MeshType1d>      DeviceType1d;

  typedef viennagrid::triangular_2d_mesh     MeshType2d;
  typedef viennashe::device<MeshType2d>      DeviceType2d;

  typedef viennagrid::tetrahedral_3d_mesh    MeshType3d;
  typedef viennashe::device<MeshType3d>      DeviceType3d;

  std::cout << "* main(): Testing mesh in 1d..." << std::endl;
  DeviceType1d device1d;
  viennashe::util::device_generation_config generator_params;
  generator_params.add_segment(0.0,  1e-6, 30);  //start at x=0, length 1e-6, 30 points
  device1d.generate_mesh(generator_params);

  if (test(device1d, "quantity_transfer_1d") != EXIT_SUCCESS)
    return EXIT_FAILURE;



  std::cout << "* main(): Testing mesh in 2d..." << std::endl;
  DeviceType2d device2d;
  device2d.load_mesh("../../examples/data/nin2d.mesh");

  if (test(device2d, "quantity_transfer_2d") != EXIT_SUCCESS)
    return EXIT_FAILURE;


  std::cout << "* main(): Testing mesh in 3d..." << std::endl;
  DeviceType3d device3d;
  device3d.load_mesh("../../examples/data/half-trigate57656.mesh");

  if (test(device3d, "quantity_transfer_3d") != EXIT_SUCCESS)
    return EXIT_FAILURE;

}

