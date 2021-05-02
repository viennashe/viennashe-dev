#ifndef VIENNASHE_UTIL_CHECKS_HPP
#define VIENNASHE_UTIL_CHECKS_HPP

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

// std
#include <iostream>
#include <limits>
#include <vector>
#include <cmath>

// viennashe
#include "viennashe/util/exception.hpp"
#include "viennashe/math/linalg_util.hpp"

#include "viennashe/log/log.hpp"
#include "viennashe/util/log_keys.h"


/** @file viennashe/util/checks.hpp
    @brief Defines a set of checker functors for micro-tests within ViennaSHE.
*/

namespace viennashe
{
  namespace util
  {


    //
    // Checker primitives (returns false if check fails)
    //

    /** @brief A checker that never complains. Can be used instead of overly picky/erroneous checkers (just for debugging purposes, of course!) */
    struct all_ok_checker
    {
      template <typename T>
      bool operator()(T const &) const { return true; }
    };


    /** @brief Checker for values inside a half-open interval [a, b), where 'a' is included, while 'b' is not. */
    template <typename ValueType = double>
    class range_checker
    {
      public:
        range_checker(ValueType const & a, ValueType const & b) : lower_bound_(a), upper_bound_(b)
        {
          if (! (a < b) )
          {
            std::stringstream ss;
            ss << "Invalid range [" << a << ", " << b << ")";
            throw invalid_range_bounds_exception(ss.str());
          }
        }

        bool operator()(ValueType const & v) const
        {
          return (v >= lower_bound_) && (v < upper_bound_);
        }

      private:
        ValueType lower_bound_;
        ValueType upper_bound_;
    };


    /** @brief Checks for NaN */
    struct nan_checker
    {
      template <typename T>
      bool operator()(T const & t) const { return t != t; }
    };

    /** @brief Checks for Inf */
    struct inf_checker
    {
      template <typename T>
      bool operator()(T const & t) const { return t - t != 0; }
    };


    /** @brief Checks if a value of type ValueType is NaN using value != value. */
    template < typename ValueType >
    bool is_NaN(const ValueType & val) { return (val != val); }

    /** @brief Checks if a value of type ValueType is Inf using (value - value) != 0. */
    template < typename ValueType >
    bool is_Inf(const ValueType & val) { return std::fabs(val - val) > 0; }

    /** @brief Checks if a value of type ValueType is negative using value < 0. */
    template < typename ValueType >
    bool is_negative(const ValueType & val) { return (val < 0); }



    /** @brief Checks a vector for valid entries (i.e. no NaN). */
    template <typename VectorType>
    void check_vector_for_valid_entries(VectorType const & vec, std::string message = "Location not specified")
    {
      for (std::size_t i=0; i<vec.size(); ++i)
      {
        if ( viennashe::util::is_NaN(vec[i]) )
        {
          std::stringstream ss;
          ss << message << ": Found NaN in vector at index " << i << std::endl;
          throw invalid_value_exception(ss.str());
        }

        if ( viennashe::util::is_Inf(vec[i]) )
        {
          std::stringstream ss;
          ss << message << ": Found infinity in vector at index " << i << std::endl;
          throw invalid_value_exception(ss.str());
        }
      }
    }


    //
    // Checkers with exceptions
    //

    /** @brief A checker class that throws a user-provided exception rather than returning false on its functor interface. */
    template <typename BasicChecker, typename ExceptionType>
    class checker_with_exception
    {
      public:
        checker_with_exception(BasicChecker const & c,
                               ExceptionType const & e) : checker_(c), exception_(e) {}

        template <typename T>
        bool operator()(T const & t) const
        {
          if (!checker_(t))
            throw exception_;

          return true;
        }

      private:
        BasicChecker const & checker_;
        ExceptionType const & exception_;
    };


    /** @brief Convenience creator routine for creating a checker with exception from checker returning a bool only */
    template <typename CheckerType, typename ExceptionType>
    checker_with_exception<CheckerType, ExceptionType> make_checker_with_exception(CheckerType const & checker, ExceptionType const & ex)
    {
      return checker_with_exception<CheckerType, ExceptionType>(checker, ex);
    }

    /** @brief Checks a matrix for being an M matrix
     *
     * To be an M-matrix, in each row the diagonal entry is the only positive entry.
     * Moreover, the modulus of the sum of the off-diagonal entries has to be smaller or equal (in at least one row strictly smaller) than the off-diagonal entry.
     *
     *  @param   A   The matrix to be checked for being an M-matrix
     */
    template <typename MatrixType>
    void m_matrix_check(MatrixType const & A)
    {
      bool is_m_matrix = true;

      for (std::size_t row_index  = 0;
                       row_index != A.size1();
                     ++row_index)
      {
        typedef typename MatrixType::const_iterator2   AlongRowIterator;
        typedef typename MatrixType::row_type          RowType;

        RowType const & row_i = A.row(row_index);

        double abs_offdiagonal = 0.0;
        double abs_diagonal = 0.0;

        for (AlongRowIterator col_it  = row_i.begin();
                              col_it != row_i.end();
                            ++col_it)
        {
          std::size_t col_index = col_it->first;
          double entry = col_it->second;

          //col_buffer.set(data_index, col_it->first);
          //elements[data_index] = col_it->second;
          //++data_index;

          if (row_index != col_index)
          {
              if (entry > 0.0)
              {
                is_m_matrix = false;
                log::warn() << "* WARNING in m_matrix_check(): Entry (" << row_index << "," << col_index << ") is positive: " << entry << std::endl;
              }
              else
                abs_offdiagonal += std::fabs(col_it->second);
          }
          else
          {
              if (entry <= 0.0)
              {
                is_m_matrix = false;
                log::warn() << "* WARNING in m_matrix_check(): Entry (" << row_index << "," << col_index << ") is non-positive: " << entry << std::endl;
              }

              abs_diagonal = entry;
          }
        }

        if (abs_offdiagonal > abs_diagonal * 1.0001)  //allow some numerical noise
        {
          is_m_matrix = false;
          log::warn() << "* WARNING in m_matrix_check(): Row " << row_index << " imbalanced! "
                    << "fabs(offdiagonals): " << abs_offdiagonal << ", diagonal: " << abs_diagonal << std::endl;
        }
      }

      if (is_m_matrix)
        log::info<log_m_matrix_check>() << "* m_matrix_check(): Matrix is M-matrix!" << std::endl;

    }





    /** @brief Checks a matrix for empty rows.
     *
     * @return -1 if everything is okay, otherwise index of first empty row
    *
    */
    template <typename MatrixType>
    long matrix_consistency_check(MatrixType const & matrix)
    {
      typedef typename MatrixType::const_iterator1       RowIterator;
      typedef typename MatrixType::const_iterator2       ColumnIterator;

      //long last_failed_index = -1; // [MB] unused variable

      for (RowIterator row_iter = matrix.begin1();
           row_iter != matrix.end1();
           ++row_iter)
      {
        bool has_row = false;
        for (ColumnIterator col_iter = row_iter.begin();
             col_iter != row_iter.end();
             ++col_iter)
        {
          if (*col_iter != 0)
            has_row = true;
        }

        if (has_row == false)
        {
          log::error() << "* FATAL ERROR in check_matrix_consistency(): System matrix has empty row " << row_iter.index1() << "!" << std::endl;
          //last_failed_index = row_iter.index1(); // [MB] row_iter.index1() has no side effects, right?
        }
      }

      return -1;
    }


    /** @brief Checks a matrix for empty rows.
     *
     * @return -1 if everything is okay, otherwise index of first empty row
    *
    */
    template <typename NumericT>
    long matrix_consistency_check(viennashe::math::sparse_matrix<NumericT> const & matrix)
    {
      typedef typename viennashe::math::sparse_matrix<NumericT>::const_iterator2   AlongRowIterator;

      for (std::size_t i = 0; i < matrix.size1(); ++i)
      {
        bool has_row = false;
        for (AlongRowIterator col_iter  = matrix.row(i).begin();
                              col_iter != matrix.row(i).end();
                            ++col_iter)
        {
          if (col_iter->second)
          {
            has_row = true;
            break;
          }
        }

        if (has_row == false)
        {
          log::error() << "* FATAL ERROR in check_matrix_consistency(): System matrix has empty row " << i << "!" << std::endl;
          //last_failed_index = row_iter.index1(); // [MB] row_iter.index1() has no side effects, right?
        }
      }

      return -1;
    }


    /** @brief Checks that the first 'num_cols' column sums of the provided
     *         matrix vanish (up to round-off) or a column is not populated with any entries at all. */
    template <typename MatrixType>
    void check_vanishing_column_sums(MatrixType const & matrix, std::size_t num_cols)
    {
      typedef typename MatrixType::const_iterator1       RowIterator;
      typedef typename MatrixType::const_iterator2       ColumnIterator;
      typedef typename MatrixType::value_type      ScalarType;

      std::size_t size = std::min<std::size_t>(matrix.size2(), num_cols);
      std::vector<ScalarType> col_sums(size);
      std::vector<ScalarType> col_maximums(size);

      //
      // Step 1: Sum into column-sum array and remember maximium in each column
      //
      for (RowIterator row_iter = matrix.begin1();
           row_iter != matrix.end1();
           ++row_iter)
      {
        if (row_iter.index1() >= size)
          break;

        for (ColumnIterator col_iter = row_iter.begin();
             col_iter != row_iter.end();
             ++col_iter)
        {
          col_sums[col_iter.index2()] += *col_iter;
          if (std::abs(*col_iter) > col_maximums[col_iter.index2()])
            col_maximums[col_iter.index2()] = std::abs(*col_iter);
        }

      }

      //
      // Step 2: Perform the tests
      //
      for (std::size_t i=0; i<col_sums.size(); ++i)
      {
        if (col_maximums[i] == 0)   //Note: If there is no inelastic scattering, column maximums are supposed to be zero!
        {
        }
        else if (std::abs(col_sums[i]) / col_maximums[i] > 1e-10)
        {
          log::error() << "* CHECK FAILED in check_vanishing_column_sums(): Column " << i << " does not vanish: " << col_sums[i] << " with maximum " << col_maximums[i] << std::endl;
          // TODO: Think about throwing an exception here? Maybe in debug mode only?
        }
      }
    }


  } //namespace util
} //namespace viennashe

#endif
