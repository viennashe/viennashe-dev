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
#include "viennagrid/viennagrid.h"


/** @brief A field containter for values of ValueType. Supports the accessor interface */
template <typename ValueType>
struct field_container
{
  typedef ValueType   value_type;

  field_container(std::size_t size, value_type default_value = value_type()) : data_(size, default_value) {}

  template <typename ElementType>
  value_type operator()(ElementType const & elem) const
  {
    return data_.at(std::size_t(viennagrid_index_from_element_id(elem)));
  }

  template <typename ElementType>
  void operator()(ElementType const & elem, value_type val)
  {
    data_.at(std::size_t(viennagrid_index_from_element_id(elem))) = val;
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

  std::vector<double> reference_field(3);
  reference_field[0] = 1.0;
  reference_field[1] = 2.0;
  reference_field[2] = 3.0;

  MeshType const & mesh = device.mesh();

  viennagrid_dimension geo_dim;
  viennagrid_mesh_geometric_dimension_get(device.mesh(), &geo_dim);

  //
  // Write normal components of field (1, 1, 0) to each edge:
  //
  std::cout << "Writing normal projections on box interfaces to facets..." << std::endl;
  viennagrid_dimension cell_dim;
  viennagrid_mesh_cell_dimension_get(mesh, &cell_dim);

  viennagrid_element_id *facets_begin, *facets_end;
  viennagrid_mesh_elements_get(mesh, cell_dim - 1, &facets_begin, &facets_end);

  field_container<double> field_component_on_edge_container(facets_end - facets_begin);

  for (viennagrid_element_id *fit  = facets_begin;
                              fit != facets_end;
                            ++fit)
  {
    viennagrid_element_id *cells_on_facet_begin, *cells_on_facet_end;
    viennagrid_element_coboundary_elements(device.mesh(), *fit, cell_dim, &cells_on_facet_begin, &cells_on_facet_end);

    // Reference orientation of facet is from first coboundary cell to second coboundary cell:

    std::vector<double> n = viennashe::util::outer_cell_normal_at_facet(device.mesh(), cells_on_facet_begin[0], *fit);
    double norm = 0;
    for (std::size_t i=0; i<n.size(); ++i)
      norm += n[i] * n[i];
    norm = std::sqrt(norm);
    for (std::size_t i=0; i<n.size(); ++i)
      n[i] /= norm; // n should be normalized already, but let's ensure this is indeed the case

    double field_in_n = reference_field[0] * n[0];
    if (geo_dim >= 2)
      field_in_n += reference_field[1] * n[1];
    if (geo_dim >= 3)
      field_in_n += reference_field[2] * n[2];
    field_component_on_edge_container(*fit, field_in_n);
  }


  //
  // Now transfer from normal projection on edges to vertices:
  //
  std::cout << "Transferring from edges to vertices..." << std::endl;
  viennagrid_element_id *cells_begin, *cells_end;
  viennagrid_mesh_elements_get(mesh, cell_dim, &cells_begin, &cells_end);

  field_container<std::vector<double> > my_cell_field_container(cells_end - cells_begin, std::vector<double>(geo_dim));

  for (viennagrid_element_id *cit  = cells_begin;
                              cit != cells_end;
                            ++cit)
  {
    /* TODO: Fix this part!
    FacetOnCellContainer facets_on_cell(*cit);

    viennashe::util::dual_box_flux_to_cell(device,
                                           *cit,                     facets_on_cell,
                                           my_cell_field_container,  field_component_on_edge_container); */

    std::vector<double> e_field(3);
    e_field[0] = my_cell_field_container(*cit)[0];
    if (geo_dim >= 2)
      e_field[1] = my_cell_field_container(*cit)[1];
    if (geo_dim >= 3)
      e_field[2] = my_cell_field_container(*cit)[2];

    //
    // Run checks:
    //
    double diff  = (e_field[0] - reference_field[0]) * (e_field[0] - reference_field[0]);
    if (geo_dim >= 2)
           diff += (e_field[1] - reference_field[1]) * (e_field[1] - reference_field[1]);
    if (geo_dim >= 3)
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

  /* TODO: Migrate to ViennaGrid 3.0
  std::cout << "Writing to VTK..." << std::endl;
  viennagrid::io::vtk_writer<MeshType> my_vtk_writer;
  my_vtk_writer.add_vector_data_on_cells(viennagrid::make_accessor<CellType>(my_cell_field_container.container()), "flux");
  my_vtk_writer(mesh, output_filename); */

  //
  // Test result:
  //

  return EXIT_SUCCESS;
}



int main()
{
  std::cout << "* main(): Testing mesh in 1d..." << std::endl;
  viennashe::device device1d;

  viennashe::util::device_generation_config generator_params;
  generator_params.add_segment(0.0,  1e-6, 30);  //start at x=0, length 1e-6, 30 points
  device1d.generate_mesh(generator_params);

  if (test(device1d, "quantity_transfer_1d") != EXIT_SUCCESS)
    return EXIT_FAILURE;



  std::cout << "* main(): Testing mesh in 2d..." << std::endl;
  viennashe::device device2d;
  device2d.load_mesh("../../examples/data/nin2d.mesh");

  if (test(device2d, "quantity_transfer_2d") != EXIT_SUCCESS)
    return EXIT_FAILURE;


  std::cout << "* main(): Testing mesh in 3d..." << std::endl;
  viennashe::device device3d;
  device3d.load_mesh("../../examples/data/half-trigate57656.mesh");

  if (test(device3d, "quantity_transfer_3d") != EXIT_SUCCESS)
    return EXIT_FAILURE;

}

