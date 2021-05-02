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

/** @example resistor1d-c.c Simple 1D simulation of a resistor using the C-library libviennashe */

#include <stdio.h>
#include <stdlib.h>

/* Main include for libviennashe */
#include "../../libviennashe/include/libviennashe.h"

/* Definition of true and false */
#define FALSE libviennashe_false
#define TRUE  libviennashe_true

/* Define the number of points to use for the resistor */
#define POINTSX 51

/** @brief Prints the names of the available quantities in a register  */
void print_quan_info(viennashe_quan_register reg)
{
  /* Initialize the variables */
  char ** names = NULL;
  viennashe_index_type num_quan = 0;
  viennashe_index_type i = 0;

  /* Get the number of available quantities */
  viennashe_get_num_cell_based(reg, &num_quan);
  printf("Number of Quanitites: %ld\n", num_quan);

  /* Retrieve the quantity names */
  names = (char**)malloc(sizeof(char*) * num_quan);
  viennashe_get_cell_based_quantity_list(reg, names);

  /* Print the names on stdout */
  printf("#############################\n");
  for (i = 0; i < num_quan; ++i)
  {
    printf("%3ld: '%s'\n", i, names[i]);
    free(names[i]);
  }
  printf("#############################\n");

  /* Do not forget: Free memory! */
  free(names); names = NULL;
}

/** @brief Prints the electrostatic potential profile to stdout  */
void print_potential(viennashe_device dev, viennashe_quan_register reg)
{
  /* Initialize the variables  */
  viennashe_index_type num_cells = 0;
  viennashe_index_type * len = NULL;
  double ** values = NULL;
  viennashe_index_type i = 0;
  viennashe_index_type j = 0;
  libviennashe_bool check = FALSE;

  /* First check if there is a quantity named "Electrostatic potential" */
  viennashe_has_cell_based_quantity(reg, "Electrostatic potential", &check);
  if(check)
  {
    /* Get the number of cells to preallocate memory */
    viennashe_get_num_cells(dev, &num_cells);

    /* Allocate memory to hold all the values */
    viennashe_prealloc_cell_based_quantity(dev, &values, &len);

    /* Retrieve the values */
    viennashe_get_cell_based_quantity(reg, "Electrostatic potential", values, len);

    /* Print the values to stdout */
    for ( i = 0; i < num_cells; ++i)
    {
      printf("%2ld => ", i);
      for ( j = 0; j < len[i]; ++j)
      {
        printf("%e", values[i][j]);
      }
      printf("\n");
    }

    /* Free previously allocated memory */
    viennashe_free_cell_based_quantity(dev, &values, &len);
  }
  else
  {
    /* If there is no electrostatic potential print an error message */
   printf("Unable to find the electrostatic potential by name ('Electrostatic potential') !\n");
  }
}


/*
 *
 */
int main()
{
  /* Prealloc needed variables */
  const double len_x    = 1e-6;
  const long   points_x = POINTSX;

  /* Material and doping */
  viennashe_material_id   matids[POINTSX-1];
  double Nd[POINTSX-1];
  double Na[POINTSX-1];
  /* Array of cell ids on which to set boundary potentials  */
  viennashe_index_type   bnd_cells[] = { 0, 49 };
  /* Array of boundary potentials */
  double bnd_pot[] = { 0.0, 1.0 };
  long   i = 0;

  /* Handles */
  viennashe_device     dev     = NULL;
  viennashe_config     conf    = NULL;
  viennashe_simulator  sim     = NULL;
  viennashe_config     conf_dd = NULL;
  viennashe_simulator  sim_dd  = NULL;
  viennashe_quan_register reg  = NULL;

  /* Init arrays */
  viennashe_get_metal_id(&matids[0]);
  viennashe_get_metal_id(&matids[points_x-2]);

  Nd[0]          = 0;
  Nd[points_x-2] = 0;
  Na[0]          = 0;
  Na[points_x-2] = 0;

  for (i = 1; i < points_x - 1; i++)
  {
    if(i < points_x - 2) viennashe_get_silicon_id(&matids[i]);
    Na[i] = 1e7;
    Nd[i] = 1e25;
  }
  /* +++++++++++ */

  /* Init ViennaSHE */
  viennashe_initalize();

  /* Create device */
  viennashe_create_1d_device(&dev, len_x, points_x);

  /* Init the device */
  viennashe_initalize_device(dev, matids, Nd, Na);
  viennashe_set_contact_potential_cells(dev, bnd_cells, bnd_pot, 2);

  /* Create the simulator configuration */
  viennashe_create_config(&conf_dd);

  /* Apply standard DD configuration */
  viennashe_config_standard_dd(conf_dd);
  viennashe_set_nonlinear_solver_config(conf_dd, viennashe_nonlinear_solver_gummel, 55, 0.5);

  /* Get the simulator instance */
  viennashe_create_simulator(&sim_dd, dev, conf_dd);

  /* This is DD, so we do not need to give any inital guess */
  /* RUN! */
  viennashe_run(sim_dd);

  viennashe_create_quantity_register(&reg, sim_dd );

  print_quan_info(reg);

  /* print_potential(dom, reg); */

  /* Write solution to CSV files for gnuplot */
  viennashe_write_to_gnuplot(reg, "Electrostatic potential",   "potential_dd.dat");
  viennashe_write_to_gnuplot(reg, "Electron density",          "n_dd.dat");
  viennashe_write_to_gnuplot(reg, "Hole density",              "p_dd.dat");

  /* Free the quantity register */
  viennashe_free_quantity_register(reg); reg = NULL;

  /*  A single SHE postprocessing step  */
  viennashe_create_config(&conf);

  /* Bipolar SHE (n and p) */
  viennashe_config_she_bipolar(conf, FALSE);
  viennashe_set_nonlinear_solver_config(conf, viennashe_nonlinear_solver_gummel, 1, 0.5);

  /* Get the simulator instance */
  viennashe_create_simulator(&sim, dev, conf);

  /* Transfer DD solution as init guess to SHE */
  viennashe_set_initial_guess_from_other_sim(sim, sim_dd);

  /* Free memory */
  viennashe_free_simulator(sim_dd); sim_dd  = NULL;
  viennashe_free_config(conf_dd);   conf_dd = NULL;

  /* RUN SHE! */
  viennashe_run(sim);

  /* Get the quantity register for the SHE-simulator */
  viennashe_create_quantity_register(&reg, sim);

  /* Print some info */
  print_quan_info(reg);

  /* Write pot, n and p to CSV files for gnuplot */
  viennashe_write_to_gnuplot(reg, "Electrostatic potential",     "potential_she.dat");
  viennashe_write_to_gnuplot(reg, "Electron density",            "n_she.dat");
  viennashe_write_to_gnuplot(reg, "Hole density",                "p_she.dat");
  /* Write the complete SHE result */
  viennashe_write_she_results_to_gnuplot(reg, viennashe_electron_id, "she_result_n.dat");
  viennashe_write_she_results_to_gnuplot(reg, viennashe_hole_id,     "she_result_p.dat");
  /* Write the average carrier energy to file */
  viennashe_write_to_gnuplot(reg, "Average electron energy SHE", "n_energy.dat");

  /* Free memory */
  viennashe_free_quantity_register(reg); reg = NULL;
  viennashe_free_simulator(sim); sim  = NULL;
  viennashe_free_config(conf);   conf = NULL;
  viennashe_free_device(dev);    dev  = NULL;

  viennashe_finalize();

  return (EXIT_SUCCESS);
}

