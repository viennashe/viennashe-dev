/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
    ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

   Additional information can be found on the project-webpage: http://viennashe.sourceforge.net

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

%inline %{ /* 1D creator for the device */
viennashe_device create_1d_device(double len_x, long points_x) {
  viennashe_device   ptr;
  viennasheErrorCode returnValue;
  returnValue = viennashe_create_1d_device(&ptr, len_x, points_x);
  if (returnValue != 0) ptr = NULL;
  return ptr;
}%}

%inline %{ /* creator for the device from file */
viennashe_device create_device_from_file(viennashe_topology_type_id topology_id, char const * filename) {
  viennashe_device   ptr;
  viennasheErrorCode returnValue;
  returnValue = viennashe_create_device_from_file(&ptr, topology_id, filename);
  if (returnValue != 0) ptr = NULL;
  return ptr;
}%}


%inline %{ /* creator for the device from arrays/lists */
viennashe_device create_device(viennashe_topology_type_id tid, PyObject * vertices, PyObject * cells, PyObject * segmentation) {
  viennashe_device   ptr = NULL;
  viennasheErrorCode returnValue = 1;
  viennashe_index_type num_vertices = 0;
  viennashe_index_type num_cells = 0;

  int dim = 0;
  int nen = 0;
  if (tid == viennashe_line_1d)          { dim = 1; nen = 2; }
  if (tid == viennashe_quadrilateral_2d) { dim = 2; nen = 4; }
  if (tid == viennashe_triangular_2d)    { dim = 2; nen = 3; }
  if (tid == viennashe_hexahedral_3d)    { dim = 3; nen = 6; }
  if (tid == viennashe_tetrahedral_3d)   { dim = 3; nen = 4; }


  double * vertices_ = NULL;
  viennashe_index_type * cells_ = NULL;
  viennashe_index_type * segmentation_ = NULL;

  if (PyList_Check(vertices))
  {
    int i = 0;
    int size = PyList_Size(vertices);
    if (size <= 0)
    {
      PyErr_SetString(PyExc_TypeError, "vertices list must contain vertices");
      return NULL;
    }
    num_vertices = size;

    vertices_ = (double *) malloc( dim * size * sizeof(double));

    for (i = 0; i < size; i++)
    {
      PyObject * co = PyList_GetItem(vertices, i);
      if (PyList_Check(co))
      {
        int j = 0;
        if (PyList_Size(co) != dim)
        {
          if(vertices_) free(vertices_);
          PyErr_SetString(PyExc_TypeError, "vertices list must contain actual coordinates. Check dimensions!");
          return NULL;
        }


        for (j = 0; j < dim; j++)
        {
          PyObject * o = PyList_GetItem(co, j);
          if (PyFloat_Check(o))
          {
            vertices_[i*dim + j] = PyFloat_AsDouble(o);
          }
          else
          {
            PyErr_SetString(PyExc_TypeError, "vertices list must contain floating point values");
            if(vertices_) free(vertices_);
            return NULL;
          }
        }
      }
      else
      {
        PyErr_SetString(PyExc_TypeError, "vertices list must contain vertices");
        if(vertices_) free(vertices_);
        return NULL;
      }
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "not a list");
    if(vertices_) free(vertices_);
    return NULL;
  }

  if (PyList_Check(cells))
  {
    int i = 0;
    int size = PyList_Size(cells);
    if (size <= 0)
    {
      PyErr_SetString(PyExc_TypeError, "cells list must contain cells");
      return NULL;
    }
    num_cells = size;

    cells_ = (viennashe_index_type * ) malloc( nen * size * sizeof(viennashe_index_type));

    for (i = 0; i < size; i++)
    {
      PyObject * co = PyList_GetItem(cells, i);
      if (PyList_Check(co))
      {
        int j = 0;
        if (PyList_Size(co) != nen)
        {
          if(vertices_) free(vertices_);
          if(cells_) free(cells_);
          PyErr_SetString(PyExc_TypeError, "cells list must contain actual vertex ids. Check dimensions!");
          return NULL;
        }

        for (j = 0; j < nen; j++)
        {
          PyObject * o = PyList_GetItem(co, j);
          if (PyInt_Check(o))
          {
            cells_[i*nen + j] = ((viennashe_index_type)PyInt_AsLong(o));
          }
          else
          {
            PyErr_SetString(PyExc_TypeError, "cells list must contain integer values");
            if(vertices_) free(vertices_);
            if(cells_) free(cells_);
            return NULL;
          }
        }
      }
      else
      {
        PyErr_SetString(PyExc_TypeError, "cells list must contain vertex indices");
        if(vertices_) free(vertices_);
        if(cells_) free(cells_);
        return NULL;
      }
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "not a list");
    if(vertices_) free(vertices_);
    if(cells_) free(cells_);
    return NULL;
  }

  if (PyList_Check(segmentation))
  {
    int i = 0;
    int size = PyList_Size(segmentation);
    if (size > 0)
    {
      segmentation_ = (viennashe_index_type * ) malloc(size * sizeof(viennashe_index_type));

      for (i = 0; i < size; i++)
      {
        PyObject * co = PyList_GetItem(segmentation, i);
        if (PyInt_Check(co))
        {
          segmentation_[i] = (viennashe_index_type)PyInt_AsLong(co);
        }
        else
        {
          PyErr_SetString(PyExc_TypeError, "segmentation list must contain segment indices");
          if(vertices_) free(vertices_);
          if(cells_) free(cells_);
          if(segmentation_) free(cells_);
          return NULL;
        }
      }
    }
  }

  returnValue = viennashe_create_device_flat(&ptr, tid, vertices_, num_vertices, cells_, num_cells, segmentation_);

  if(vertices_) free(vertices_);
  if(cells_) free(cells_);

  if (returnValue != 0) { ptr = NULL; }
  return ptr;
}%}


%newobject create_1d_device;
%newobject create_device_from_file;
%newobject create_device;
%delobject viennashe_free_device;




