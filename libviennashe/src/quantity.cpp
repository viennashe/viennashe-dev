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
#include "libviennashe/src/quantity_register.hpp"

// C includes
#include "libviennashe/include/quantity.h"

#include <cstring>


#ifdef	__cplusplus
extern "C"
{
#endif

viennasheErrorCode viennashe_create_quantity_register(viennashe_quan_register * reg, viennashe_simulator sim)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(sim,2,"sim");

    // Get simulator
    viennashe_simulator_impl * int_sim = sim;
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! viennashe_create_quantity_register(): The simulator (sim) must be valid!" << std::endl;
      return 2;
    }

    libviennashe::quan_register_internal * int_reg = new libviennashe::quan_register_internal;
    int_reg->int_sim = int_sim; // Make sure to register the quans for a certain simulator!


    //
    // Register Quantities
    if (int_sim->stype == libviennashe::meshtype::line_1d)
      libviennashe::register_quans(*(int_reg->int_sim->sim1d),  *int_reg);
    else if (int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
      libviennashe::register_quans(*(int_reg->int_sim->simq2d), *int_reg);
    else if (int_sim->stype == libviennashe::meshtype::triangular_2d)
      libviennashe::register_quans(*(int_reg->int_sim->simt2d), *int_reg);
    else if (int_sim->stype == libviennashe::meshtype::hexahedral_3d)
      libviennashe::register_quans(*(int_reg->int_sim->simh3d), *int_reg);
    else if (int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
      libviennashe::register_quans(*(int_reg->int_sim->simt3d), *int_reg);
    else
    {
      viennashe::log::error() << "ERROR! viennashe_create_quantity_register(): Malconfigured simulator!" << std::endl;
      return 2;
    }

    *reg = reinterpret_cast<viennashe_quan_register>(int_reg);
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_create_quantity_register(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_free_quantity_register(viennashe_quan_register reg)
{
  try
  {
    if (reg != NULL)
    {
      delete reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    }
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_free_quantity_register(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


viennasheErrorCode viennashe_get_num_cell_based(viennashe_quan_register reg, viennashe_index_type * num)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(num,2,"num");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    *num = static_cast<viennashe_index_type>(int_reg->cell_based.size());
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_num_vertex_based(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


viennasheErrorCode viennashe_get_cell_based_quantity_list(viennashe_quan_register reg, char ** names)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(names,2,"names");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    std::vector<std::string> vnames = int_reg->cell_based.get_names();
    for (std::size_t i = 0; i < vnames.size(); ++i)
    {
      // No "new" statement here ... since valgrind would complain ...
      names[i] = (char*) malloc(sizeof(char) * (vnames[i].size()+1));
      std::strncpy(names[i], vnames[i].c_str(), vnames[i].size()+1);
    }
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_vertex_based_quantity_list(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

VIENNASHE_EXPORT viennasheErrorCode viennashe_has_cell_based_quantity(viennashe_quan_register reg, char const * name, libviennashe_bool * exists)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(name,2,"name");
    CHECK_ARGUMENT_FOR_NULL(exists,3,"exists");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);

    // Translate bool to LIBVIENNASHE_BOOL "by hand" just to be sure
    if(int_reg->cell_based.has_quan(std::string(name)) == true)
      *exists = libviennashe_true;
    else
      *exists = libviennashe_false;

  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_vertex_based_quantity(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

viennasheErrorCode viennashe_get_cell_based_quantity(viennashe_quan_register reg, char const * name, double ** values, viennashe_index_type * len)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(name,2,"name");
    CHECK_ARGUMENT_FOR_NULL(values,3,"values");
    CHECK_ARGUMENT_FOR_NULL(len,4,"len");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);

    if (!int_reg->cell_based.has_quan(std::string(name)))
    {
      viennashe::log::error() << "ERROR: viennashe_get_vertex_based_quantity(): The quantity '" << name << "' does not exist." << std::endl;
      return 2;
    }
    // fill value array
    int_reg->cell_based.get(std::string(name)).fill(values, len);

  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_vertex_based_quantity(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_edf(viennashe_quan_register reg, viennashe_carrier_ids ctype,
                                                          double ** energies, double ** values, viennashe_index_type * len)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(len,5,"len");
    CHECK_ARGUMENT_FOR_NULL(energies,3,"energies");
    CHECK_ARGUMENT_FOR_NULL(values,4,"values");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    // Get simulator
    viennashe_simulator_impl * int_sim = int_reg->int_sim;
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_edf(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    viennashe::carrier_type_id arg_ctype = (ctype == viennashe_electron_id) ? viennashe::ELECTRON_TYPE_ID : viennashe::HOLE_TYPE_ID;

    if (int_sim->stype == libviennashe::meshtype::line_1d)
      libviennashe::she_fill_edf(*(int_reg->int_sim->sim1d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
      libviennashe::she_fill_edf(*(int_reg->int_sim->simq2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::triangular_2d)
      libviennashe::she_fill_edf(*(int_reg->int_sim->simt2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::hexahedral_3d)
      libviennashe::she_fill_edf(*(int_reg->int_sim->simh3d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
      libviennashe::she_fill_edf(*(int_reg->int_sim->simt3d), arg_ctype, energies, values, len);
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_edf(): Malconfigured simulator!" << std::endl;
      return 2;
    }


  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_she_edf(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_dos(viennashe_quan_register reg, viennashe_carrier_ids ctype,
                                                          double ** energies, double ** values, viennashe_index_type * len)
{
  try
  {
    if (reg == NULL)
    {
      viennashe::log::error() << "ERROR: viennashe_get_she_dos(): The register pointer (reg) must not be NULL." << std::endl;
      return 1;
    }
    if (len == NULL) { return 5; }
    if (energies == NULL) { return 3; }
    if (values == NULL) { return 4; }
    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    // Get simulator
    viennashe_simulator_impl * int_sim = int_reg->int_sim;
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_dos(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    viennashe::carrier_type_id arg_ctype = (ctype == viennashe_electron_id) ? viennashe::ELECTRON_TYPE_ID : viennashe::HOLE_TYPE_ID;

    if (int_sim->stype == libviennashe::meshtype::line_1d)
      libviennashe::she_fill_dos(*(int_reg->int_sim->sim1d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
      libviennashe::she_fill_dos(*(int_reg->int_sim->simq2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::triangular_2d)
      libviennashe::she_fill_dos(*(int_reg->int_sim->simt2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::hexahedral_3d)
      libviennashe::she_fill_dos(*(int_reg->int_sim->simh3d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
      libviennashe::she_fill_dos(*(int_reg->int_sim->simt3d), arg_ctype, energies, values, len);
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_dos(): Malconfigured simulator!" << std::endl;
      return 2;
    }


  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_she_dos(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}

VIENNASHE_EXPORT viennasheErrorCode viennashe_get_she_group_velocity(viennashe_quan_register reg, viennashe_carrier_ids ctype,
                                                                     double ** energies, double ** values, viennashe_index_type * len)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(reg,1,"reg");
    CHECK_ARGUMENT_FOR_NULL(len,5,"len");
    CHECK_ARGUMENT_FOR_NULL(energies,3,"energies");
    CHECK_ARGUMENT_FOR_NULL(values,4,"values");

    libviennashe::quan_register_internal * int_reg = reinterpret_cast<libviennashe::quan_register_internal *>(reg);
    // Get simulator
    viennashe_simulator_impl * int_sim = int_reg->int_sim;
    if (!int_sim->is_valid())
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_group_velocity(): The simulator (sim) must be valid!" << std::endl;
      return 1;
    }

    viennashe::carrier_type_id arg_ctype = (ctype == viennashe_electron_id) ? viennashe::ELECTRON_TYPE_ID : viennashe::HOLE_TYPE_ID;

    if (int_sim->stype == libviennashe::meshtype::line_1d)
      libviennashe::she_fill_group_velocity(*(int_reg->int_sim->sim1d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::quadrilateral_2d)
      libviennashe::she_fill_group_velocity(*(int_reg->int_sim->simq2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::triangular_2d)
      libviennashe::she_fill_group_velocity(*(int_reg->int_sim->simt2d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::hexahedral_3d)
      libviennashe::she_fill_group_velocity(*(int_reg->int_sim->simh3d), arg_ctype, energies, values, len);
    else if (int_sim->stype == libviennashe::meshtype::tetrahedral_3d)
      libviennashe::she_fill_group_velocity(*(int_reg->int_sim->simt3d), arg_ctype, energies, values, len);
    else
    {
      viennashe::log::error() << "ERROR! viennashe_get_she_group_velocity(): Malconfigured simulator!" << std::endl;
      return 2;
    }


  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_get_she_group_velocity(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_prealloc_cell_based_quantity(viennashe_device dev, double *** uarray, viennashe_index_type ** len)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");

    viennashe_index_type num_vertices = 0;
    if (viennashe_get_num_vertices(dev, &num_vertices) != 0 || num_vertices <= 0)
    {
      viennashe::log::error() << "ERROR: viennashe_prealloc_vertex_based_quantity(): The mesh is malconfigured." << std::endl;
      return 1;
    }

    if ( *uarray == NULL )
    {
      // allocate
      *uarray = new double*[num_vertices];
      // NULLify
      for (size_t i = 0; i < num_vertices; ++i)
      {
        (*uarray)[i] = NULL;
      } // for each vertex
    }
    else
    {
      viennashe::log::error() << "ERROR: viennashe_prealloc_vertex_based_quantity():  The user array (uarray) is already preallocated (uarray != NULL)." << std::endl;
      return 2;
    }

    if ( *len == NULL )
    {
      // allocate
      *len = new viennashe_index_type[num_vertices];
      // NULLify
      for (size_t i = 0; i < num_vertices; ++i)
      {
        (*len)[i] = 0;
      } // for each vertex
    }
    else
    {
      viennashe::log::error() << "ERROR: viennashe_prealloc_vertex_based_quantity():  The user array (len) is already preallocated (len != NULL)." << std::endl;
      return 3;
    }


  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_prealloc_vertex_based_quantity(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


VIENNASHE_EXPORT viennasheErrorCode viennashe_free_cell_based_quantity(viennashe_device dev, double *** uarray, viennashe_index_type ** len)
{
  try
  {
    CHECK_ARGUMENT_FOR_NULL(dev,1,"dev");

    viennashe_index_type num_vertices = 0;
    if (viennashe_get_num_vertices(dev, &num_vertices) != 0 || num_vertices <= 0)
    {
      viennashe::log::error() << "ERROR: viennashe_free_vertex_based_quantity(): The mesh is malconfigured." << std::endl;
      return 1;
    }

    if ((*uarray) != NULL)
    {
      for (size_t i = 0; i < num_vertices; ++i)
      {
        if ((*uarray)[i] != NULL)
        {
          delete (*uarray)[i];
          (*uarray)[i] = NULL;
        }
      } // for each vertex
      delete *uarray;
      *uarray = NULL;
    }

    if ((*len) != NULL)
    {
      delete *len;
      *len = NULL;
    }
  }
  catch (...)
  {
    viennashe::log::error() << "ERROR: viennashe_free_vertex_based_quantity(): UNKOWN ERROR!" << std::endl;
    return -1;
  }
  return 0;
}


#ifdef	__cplusplus
}
#endif


