/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
    ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

   Additional information can be found on the project-webpage: http://viennashe.sourceforge.net

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

%module pyviennashe
%include "typemaps.i"
%{
/* Includes the header in the wrapper code */
#include "libviennashe.h"
%}

/*************************************************/
/*              RENAME SECTION                   */

// Rule for output of boolean values
%apply int * OUTPUT { libviennashe_bool * enabled };

%rename("%(strip:[viennashe_])s") ""; // Remove the viennashe_ prefix everywhere !

// Rule for all material id getters
%apply viennashe_material_id * OUTPUT { viennashe_material_id * id };

// Rule for all mesh info getters
%apply viennashe_index_type * OUTPUT { viennashe_index_type * num };

// Rule for dimension output
%apply viennashe_index_type * OUTPUT { viennashe_index_type * dim };

// Rule for solver type output
%apply int * OUTPUT { viennashe_nonlinear_solver_id * sol_id }; // note: enum values are integers ...

// Rule for max iters output
%apply long * OUTPUT { long * max_iters };

// Rule for damping factor output
%apply double * OUTPUT { double * damping };


// We have to wrap this manually ...
%ignore viennashe_get_grid;

/* Special ignore rules for device creation and destruction */
%ignore viennashe_create_device;
%ignore viennashe_create_device_flat;
%ignore viennashe_create_device_from_file;
%ignore viennashe_create_1d_device;
%ignore viennashe_free_device;

// CTOR rules for the config
%ignore viennashe_create_config;
%ignore viennashe_free_config;

// CTOR rules for the simulator
%ignore viennashe_create_simulator;
%ignore viennashe_free_simulator;

// CTOR rules for the quantity register
%ignore viennashe_create_quantity_register;
%ignore viennashe_free_quantity_register;


/*************************************************/

// This maps the python boolean to libviennashe_bool
%typemap(in)  libviennashe_bool { if($input != (PyObject *)Py_True) {  $1 = libviennashe_false; } else { $1 = libviennashe_true;} }


// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if input is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
        PyErr_SetString(PyExc_TypeError, "list must contain strings");
        free($1);
        return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError, "not a list");
    return NULL;
  }
}
// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** { free((char *) $1); }


%include "devicetypemaps.i"


/**********************************************/
/* Parse the header file to generate wrappers */
%include "libviennashe.h"

/**********************************************/


/*********************************/
/* Wrap creators and destructors */


// Device Wrappers
%include "device.i"


%inline %{ /* creator for the config */
viennashe_config create_config() {
  viennashe_config   ptr;
  viennasheErrorCode returnValue;
  returnValue = viennashe_create_config(&ptr);
  if (returnValue != 0) ptr = NULL;
  return ptr;
}%}
%newobject create_config;
%delobject viennashe_free_config;



%inline %{ /* creator for the simulator */
viennashe_simulator create_simulator(viennashe_device dev, viennashe_config conf) {
  viennashe_simulator   ptr;
  viennasheErrorCode returnValue;
  returnValue = viennashe_create_simulator(&ptr, dev, conf);
  if (returnValue != 0) ptr = NULL;
  return ptr;
}%}
%newobject create_simulator;
%delobject viennashe_free_simulator;


%inline %{ /* creator for the quantity register */
viennashe_quan_register create_quantity_register(viennashe_simulator sim) {
  viennashe_quan_register   ptr;
  viennasheErrorCode returnValue;
  returnValue = viennashe_create_quantity_register(&ptr, sim);
  if (returnValue != 0) ptr = NULL;
  return ptr;
}%}
%newobject create_quantity_register;
%delobject viennashe_free_quantity_register;




/*********************************/




