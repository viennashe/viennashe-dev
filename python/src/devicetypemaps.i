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

%typemap(in) (viennashe_material_id * material_ids) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int size = PyList_Size($input);
  $1 = (viennashe_material_id *) malloc(size * sizeof(viennashe_material_id));
  for (i = 0; i < size; i++)
  {
    PyObject *s = PyList_GetItem($input, i);
    if (!PyInt_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integers");
        return NULL;
    }
    $1[i] = PyInt_AsUnsignedLongMask(s);
  }
}
%typemap(freearg) (viennashe_material_id * material_ids) { if ($1) free($1); }


%typemap(in) (double * doping_n) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int size = PyList_Size($input);
  $1 = (double *) malloc(size * sizeof(double));
  for (i = 0; i < size; i++)
  {
    PyObject *s = PyList_GetItem($input, i);
    if (!PyFloat_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be floats");
        return NULL;
    }
    $1[i] = PyFloat_AsDouble(s);
  }
}
%typemap(freearg) (double * doping_n) { if ($1) free($1); }

%typemap(in) (double * doping_p) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int size = PyList_Size($input);
  $1 = (double *) malloc(size * sizeof(double));
  for (i = 0; i < size; i++)
  {
    PyObject *s = PyList_GetItem($input, i);
    if (!PyFloat_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be floats");
        return NULL;
    }
    $1[i] = PyFloat_AsDouble(s);
  }
}
%typemap(freearg) (double * doping_p) { if ($1) free($1); }

%typemap(in) (viennashe_index_type * cell_ids) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int size = PyList_Size($input);
  $1 = (viennashe_index_type *) malloc(size * sizeof(viennashe_index_type));
  for (i = 0; i < size; i++)
  {
    PyObject *s = PyList_GetItem($input, i);
    if (!PyInt_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integers");
        return NULL;
    }
    $1[i] = PyInt_AsUnsignedLongMask(s);
  }
}
%typemap(freearg) (viennashe_index_type * cell_ids) { if ($1) free($1); }


%typemap(in) (double * values) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int size = PyList_Size($input);
  $1 = (double *) malloc(size * sizeof(double));
  for (i = 0; i < size; i++)
  {
    PyObject *s = PyList_GetItem($input, i);
    if (!PyFloat_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be floats");
        return NULL;
    }
    $1[i] = PyFloat_AsDouble(s);
  }
}
%typemap(freearg) (double * values) { if ($1) free($1); }


%apply double * OUTPUT { double * x };
%apply double * OUTPUT { double * y };
%apply double * OUTPUT { double * z };




%typemap(in,numinputs=0) viennashe_index_type * vertex_id_list (viennashe_index_type * list)
{
  if (!PyArg_ParseTuple(args,(char *)"OO:get_nth_cell",&obj0,&obj1)) SWIG_fail;
  res1 = SWIG_ConvertPtr(obj0, &argp1,SWIGTYPE_p_viennashe_device_impl, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "get_nth_cell" "', argument " "1"" of type '" "viennashe_device""'");
  }
  arg1 = (viennashe_device)(argp1);

  viennashe_index_type nen = 0;
  viennashe_get_num_vertices_per_cell(arg1, &nen);
  $1 = malloc(sizeof(viennashe_index_type)*nen);
}
%typemap(argout) viennashe_index_type * vertex_id_list {
  int i = 0;
  viennashe_index_type nen = 0;
  viennashe_get_num_vertices_per_cell(arg1, &nen);

  PyObject * o = PyList_New(nen);
  for (i = 0; i < nen; i++)
  {
    PyList_SetItem(o, i, PyInt_FromLong($1[i]));
  }
  $result = SWIG_Python_AppendOutput($result, o);

  if($1) free($1);
}


