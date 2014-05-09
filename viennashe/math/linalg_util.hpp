#ifndef VIENNASHE_MATH_LINALG_UTIL_HPP
#define VIENNASHE_MATH_LINALG_UTIL_HPP

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


///// Various utility functions for systems of linear equations ///////////////

#include <vector>
#include <map>
#include <cmath>
#include <assert.h>

#include "viennashe/forwards.h"

/** @file viennashe/math/linalg_util.hpp
    @brief Implementation of various utilities related to linear algebra
*/

namespace viennashe
{
  namespace math
  {

    template <typename VectorType>
    typename VectorType::value_type norm_2(VectorType const & v)
    {
      typename VectorType::value_type result = 0;
      for (std::size_t i=0; i<v.size(); ++i)
        result += v[i] * v[i];
      return std::sqrt(result);
    }

    /** @brief Helper routine for the generic vector subtraction z = x - y;
      *
      * Does not use any expression template magic. This, however, is okay, as this function is in no way performance critical.
      */
    template <typename VectorType>
    VectorType subtract(VectorType const & x, VectorType const & y)
    {
      assert(x.size() == y.size() && bool("Size mismatch for subtracting y from x"));

      VectorType result(x.size());
      for (std::size_t i=0; i<x.size(); ++i)
        result[i] = x[i] - y[i];
      return result;
    }

    /** @brief Sparse matrix class based on a vector of binary trees for holding the entries.
      *
      * Introduced to overcome the scalability issues of ublas::compressed_matrix
      */
    template <typename NumericT>
    class sparse_matrix
    {
      typedef sparse_matrix<NumericT>  self_type;

    public:
      typedef std::size_t    size_type;
      typedef size_type      SizeType;

    private:
      typedef std::map<size_type, NumericT>     MapType;
      typedef std::vector<MapType>              DataContainerType;

    public:
      typedef typename MapType::iterator                      iterator2;
      typedef typename MapType::const_iterator          const_iterator2;

      typedef MapType     row_type;
      typedef MapType     RowType;

      sparse_matrix() : zero_(0) {}
      sparse_matrix(size_type num_rows, size_type num_cols) : data_(num_rows), zero_(0)
      {
        assert(num_cols == num_rows && bool("Only square matrices supported."));
        (void)num_cols;
      }

      size_type size1() const { return data_.size(); }
      size_type size2() const { return data_.size(); }

      /** @brief Non-const access to the entries of the matrix. Use this for assembly. */
      NumericT & operator()(size_type i, size_type j)
      {
        assert(i < data_.size() && j < data_.size() && bool("Access to matrix out of range!"));
        return data_.at(i)[j];
      }

      NumericT const & operator()(size_type i, size_type j) const
      {
        assert(i < data_.size() && j < data_.size() && bool("Access to matrix out of range!"));
        const_iterator2 it = data_.at(i).find(j);
        if (it == data_.at(i).end())
          return zero_;
        else
          return it->second;
      }

      MapType const & row(size_t i) const { return data_.at(i); }
      MapType       & row(size_t i)       { return data_.at(i);   }

      size_type nnz() const
      {
        size_type result = 0;
        for (size_type i=0; i<data_.size(); ++i)
          result += data_[i].size();
        return result;
      }

      self_type trans() const
      {
        self_type A_trans(data_.size(), data_.size());

        for (size_type i=0; i<data_.size(); ++i)
        {
          for (const_iterator2 it = data_[i].begin(); it != data_[i].end(); ++it)
            A_trans(it->first, i) = it->second;
        }
        return A_trans;
      }

      self_type operator*(NumericT factor) const
      {
        self_type result = *this;
        for (size_type i=0; i<data_.size(); ++i)
        {
          for (const_iterator2 it = data_[i].begin(); it != data_[i].end(); ++it)
            result(it->first, i) *= factor;
        }
        return result;
      }

      self_type & operator+=(self_type B)
      {
        for (size_type i=0; i<B.size1(); ++i)
        {
          row_type row_i = B.row(i);
          for (const_iterator2 it = row_i.begin(); it != row_i.end(); ++it)
            data_[i][it->first] += it->second;
        }
        return *this;
      }

    private:
      DataContainerType  data_;
      NumericT           zero_;   //helper member for implementing operator() const
    };


    /** @brief Normalizes an equation system such that all diagonal entries are non-negative, and such that all 2-norms of the rows are unity. */
    template <typename NumericT, typename VectorT>
    VectorT row_normalize_system(sparse_matrix<NumericT> & system_matrix, VectorT & rhs)
    {
      typedef typename sparse_matrix<NumericT>::row_type     RowType;
      typedef typename sparse_matrix<NumericT>::iterator2    AlongRowIterator;

      VectorT scale_factors(rhs.size());

      for (std::size_t i=0; i<system_matrix.size1(); ++i)
      {
        RowType & row_i = system_matrix.row(i);

        // obtain current norm of row:
        double row_norm = 0.0;
        for (AlongRowIterator iter  = row_i.begin();
                              iter != row_i.end();
                            ++iter)
        {
          row_norm += iter->second * iter->second;
        }

        row_norm = sqrt(row_norm);

        // normalize such that diagonal entry becomes positive:
        if (system_matrix(i, i) < 0.0)
          row_norm *= -1.0;

        // scale row accordingly:
        for (AlongRowIterator iter  = row_i.begin();
                              iter != row_i.end();
                            ++iter)
        {
          iter->second /= row_norm;
        }
        rhs[i] /= row_norm;

        // remember scaling factor:
        scale_factors[i] = 1.0 / row_norm;
      }

      return scale_factors;
    }

    /** @brief Converts a sparse matrix (ublas, ViennaCL) to a string. Handy debugging facility. */
    template <typename NumericT>
    std::ostream & operator<<(std::ostream & os, sparse_matrix<NumericT> const & system_matrix)
    {
      typedef typename sparse_matrix<NumericT>::const_iterator2     AlongRowIterator;

      for (std::size_t i=0; i<system_matrix.size1(); ++i)
      {
        os << std::endl << "Row " << i << ": ";
        for (AlongRowIterator iter  = system_matrix.row(i).begin();
                              iter != system_matrix.row(i).end();
                            ++iter)
        {
          os << "(" << iter->first << ", " << iter->second << "),  ";
        }
      }

      return os;
    }


    /** @brief Computes A * x  for a sparse A and a vector x.
      *
      * Does not apply expression template fancyness, but this is not a performance-critical routine anyway...
     */
    template <typename NumericT, typename VectorType>
    VectorType prod(sparse_matrix<NumericT> const & system_matrix,
                    VectorType const & x)
    {
      typedef typename sparse_matrix<NumericT>::const_iterator2    Iterator2;
      typedef typename sparse_matrix<NumericT>::row_type           RowType;

      VectorType result(system_matrix.size1());
      for (std::size_t i=0; i<system_matrix.size1(); ++i)
      {
        RowType const & row_i = system_matrix.row(i);
        double val = 0.0;
        for (Iterator2 iter2  = row_i.begin();
                       iter2 != row_i.end();
                     ++iter2)
        {
          val += iter2->second * x[iter2->first];
        }
        result[i] = val;
      }

      return result;
    }



    //////////////// Dense stuff //////////////////

    template <typename NumericT>
    class dense_matrix
    {
      typedef dense_matrix<NumericT>     self_type;
    public:
      typedef std::size_t       size_type;
      typedef size_type         SizeType;

      dense_matrix(size_type rows, size_type cols) : data_(rows * cols + 1), size1_(rows), size2_(cols)
      {
        assert(rows == cols && bool("Only square dense matrices supported!"));
      }

      size_type size1() const { return size1_; }
      size_type size2() const { return size2_; }

      NumericT const & operator()(size_type i, size_type j) const
      {
        assert(i < size1_ && j < size2_ && bool("Index out of bounds!"));
        return data_.at(i * size2_ + j);
      }

      NumericT & operator()(size_type i, size_type j)
      {
        assert(i < size1_ && j < size2_ && bool("Index out of bounds!"));
        return data_.at(i * size2_ + j);
      }

    private:
      std::vector<NumericT>   data_;
      size_type               size1_;
      size_type               size2_;
    };

  } //namespace math
}//namespace viennashe
#endif
