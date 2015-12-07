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

#include <iostream>
#include <cstdlib>
#include <map>
#include <vector>
#include <string>

#include "viennashe/forwards.h"

#include "viennashe/boxmesh/boxmesh.hpp"

struct viennashe_box_quantity_impl;


struct viennashe_boxmesh_impl
{
  typedef std::map<std::string, viennashe_box_quantity_impl>   quantity_map_type;

  viennashe_boxmesh_impl() : ref_counter_(1) {}

  std::size_t  ref_counter_;

  std::vector<viennashe_boxmesh_element_id>               boxes_;
  std::vector<std::vector<viennashe_boxmesh_element_id> > interfaces_on_boxes_;

  std::vector<viennashe_boxmesh_element_id>                                           interfaces_;
  std::vector<std::pair<viennashe_boxmesh_element_id, viennashe_boxmesh_element_id> > boxes_on_interfaces_;

  quantity_map_type         box_quantities_;
  quantity_map_type   interface_quantities_;
};

viennashe_error viennashe_boxmesh_create(viennashe_boxmesh *boxmesh)
{
  *boxmesh = new viennashe_boxmesh_impl();

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_retain(viennashe_boxmesh boxmesh)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  boxmesh->ref_counter_ += 1;

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_release(viennashe_boxmesh boxmesh)
{
  if (boxmesh->ref_counter_ > 1)
    boxmesh->ref_counter_ -= 1;
  else if (boxmesh->ref_counter_ == 1)
    delete boxmesh;
  else
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  return VIENNASHE_BOXMESH_SUCCESS;
}




viennashe_error viennashe_boxmesh_box_add(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id *box_id)
{
  viennashe_boxmesh_element_id new_id = static_cast<viennashe_boxmesh_element_id>(boxmesh->boxes_.size());
  boxmesh->boxes_.push_back(*box_id);

  if (box_id)
    *box_id = new_id;

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_interface_add(viennashe_boxmesh boxmesh,
                                                viennashe_boxmesh_element_id first_box,
                                                viennashe_boxmesh_element_id second_box,
                                                viennashe_boxmesh_element_id *interface_id)
{
  if (first_box == second_box)
    return VIENNASHE_BOXMESH_ERROR_INVALID_BOX_OR_INTERFACE;

  viennashe_boxmesh_element_id new_id = static_cast<viennashe_boxmesh_element_id>(boxmesh->interfaces_.size());
  boxmesh->interfaces_.push_back(new_id);
  boxmesh->boxes_on_interfaces_.push_back(std::make_pair(first_box, second_box));
  boxmesh->interfaces_on_boxes_[first_box].push_back(new_id);
  boxmesh->interfaces_on_boxes_[second_box].push_back(new_id);

  if (interface_id)
    *interface_id = new_id;

  return VIENNASHE_BOXMESH_SUCCESS;
}


viennashe_error viennashe_boxmesh_boxes_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id **boxes_begin, viennashe_boxmesh_element_id **boxes_end)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  *boxes_begin = &(boxmesh->boxes_[0]);
  *boxes_end   = *boxes_begin + boxmesh->boxes_.size();

  return VIENNASHE_BOXMESH_SUCCESS;
}


viennashe_error viennashe_boxmesh_interfaces_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id **interfaces_begin, viennashe_boxmesh_element_id **interfaces_end)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  *interfaces_begin = &(boxmesh->interfaces_[0]);
  *interfaces_end   = *interfaces_begin + boxmesh->interfaces_.size();

  return VIENNASHE_BOXMESH_SUCCESS;
}


viennashe_error viennashe_boxmesh_box_interfaces_get(viennashe_boxmesh boxmesh,
                                                     viennashe_boxmesh_element_id box,
                                                     viennashe_boxmesh_element_id **interfaces_begin,
                                                     viennashe_boxmesh_element_id **interfaces_end)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  if (box < 0 || std::size_t(box) >= boxmesh->boxes_.size())
    return VIENNASHE_BOXMESH_ERROR_INVALID_BOX_OR_INTERFACE;

  *interfaces_begin = &(boxmesh->interfaces_on_boxes_[box][0]);
  *interfaces_end   = *interfaces_begin + boxmesh->interfaces_on_boxes_[box].size();

  return VIENNASHE_BOXMESH_SUCCESS;
}


viennashe_error viennashe_boxmesh_interface_boxes_get(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_id interface, viennashe_boxmesh_element_id *first_box, viennashe_boxmesh_element_id *second_box)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  if (interface < 0 || std::size_t(interface) >= boxmesh->interfaces_.size())
    return VIENNASHE_BOXMESH_ERROR_INVALID_BOX_OR_INTERFACE;

  *first_box    = boxmesh->boxes_on_interfaces_[interface].first;
  *second_box   = boxmesh->boxes_on_interfaces_[interface].second;

  return VIENNASHE_BOXMESH_SUCCESS;
}


////////////////// Quantity ///////////////////


struct viennashe_box_quantity_impl
{
  viennashe_box_quantity_impl() : ref_counter_(1), data_(NULL) {}

  ~viennashe_box_quantity_impl() { if (data_) free(data_); }

  std::size_t                     ref_counter_;
  viennashe_boxmesh_element_type  where_defined_;
  viennashe_boxmesh               boxmesh_;
  std::string                     name_;
  std::size_t                     bytes_per_datum_;
  void *data_;
};


viennashe_error viennashe_boxmesh_quantity_create(viennashe_boxmesh boxmesh, viennashe_boxmesh_quantity *quantity, viennashe_boxmesh_element_type where_defined, const char *name, size_t bytes_per_datum)
{
  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  viennashe_error err = viennashe_boxmesh_quantity_get_by_name(boxmesh, where_defined, name, quantity);
  if (err)
    return err;

  if (quantity)
    return VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DEFINED;

  *quantity = new viennashe_box_quantity_impl();

  if (where_defined == VIENNASHE_BOXMESH_BOX)
    (*quantity)->data_ = malloc(bytes_per_datum * boxmesh->boxes_.size());
  else if (where_defined == VIENNASHE_BOXMESH_INTERFACE)
    (*quantity)->data_ = malloc(bytes_per_datum * boxmesh->interfaces_.size());
  else
    return VIENNASHE_BOXMESH_ERROR_INVALID_ELEMENT_TYPE;

  (*quantity)->where_defined_   = where_defined;
  (*quantity)->boxmesh_         = boxmesh;
  (*quantity)->name_            = name;
  (*quantity)->bytes_per_datum_ = bytes_per_datum;

  return VIENNASHE_BOXMESH_SUCCESS;
}


viennashe_error viennashe_boxmesh_quantity_name_get(viennashe_boxmesh_quantity quantity, const char **name)
{
  if (quantity->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DESTROYED;

  *name = quantity->name_.c_str();

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_quantity_datum_size_get(viennashe_boxmesh_quantity quantity, size_t *bytes_per_datum)
{
  if (quantity->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DESTROYED;

  *bytes_per_datum = quantity->bytes_per_datum_;

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_quantity_datum_get(viennashe_boxmesh_quantity quantity, viennashe_boxmesh_element_id box_or_interface, void **datum)
{
  if (quantity->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_QUANTITY_ALREADY_DESTROYED;

  *datum = reinterpret_cast<char*>(quantity->data_) + box_or_interface * quantity->bytes_per_datum_;

  return VIENNASHE_BOXMESH_SUCCESS;
}

viennashe_error viennashe_boxmesh_quantity_get_by_name(viennashe_boxmesh boxmesh, viennashe_boxmesh_element_type assoc, const char *name, viennashe_boxmesh_quantity *quantity)
{
  typedef typename viennashe_boxmesh_impl::quantity_map_type::iterator   ConstIterator;

  if (boxmesh->ref_counter_ < 1)
    return VIENNASHE_BOXMESH_ERROR_BOXMESH_ALREADY_DESTROYED;

  quantity = NULL;

  ConstIterator it;
  if (assoc == VIENNASHE_BOXMESH_BOX)
  {
    it = boxmesh->box_quantities_.find(name);
    if (it != boxmesh->box_quantities_.end())
      *quantity = &(it->second);
  }
  else if (assoc == VIENNASHE_BOXMESH_INTERFACE)
  {
    it = boxmesh->interface_quantities_.find(name);
    if (it != boxmesh->interface_quantities_.end())
      *quantity = &(it->second);
  }
  else
    return VIENNASHE_BOXMESH_ERROR_INVALID_ELEMENT_TYPE;

  return VIENNASHE_BOXMESH_SUCCESS;
}


