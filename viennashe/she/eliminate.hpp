#ifndef VIENNASHE_SHE_ELIMINATE_HPP
#define VIENNASHE_SHE_ELIMINATE_HPP

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
#include <cassert>

// viennashe
#include "viennashe/she/exception.hpp"
#include "viennashe/log/log.hpp"
#include "viennashe/she/log_keys.h"

/** @file viennashe/she/eliminate.hpp
    @brief Implements the elimination and the recovery of odd-order unknowns from the discrete system of equations
 */

namespace viennashe
{
  namespace she
  {

    namespace detail
    {

      template <typename FullMatrixType, typename VectorType>
      void combine_lines(FullMatrixType & system_matrix, VectorType & rhs,
                         std::size_t line1, std::size_t line_to_add, std::size_t column_to_zero)
      {

        typedef typename FullMatrixType::row_type    RowType;
        typedef typename FullMatrixType::iterator2   AlongRowIterator;

        const double s = system_matrix(line1, column_to_zero) / system_matrix(line_to_add, line_to_add);

        if (!s) return;

        RowType & values_to_add = system_matrix.row(line_to_add); // 'line_to_add'

        // multiply 'line_to_add' with 's' and subtract it from the current line
        rhs[line1] = rhs[line1] - rhs[line_to_add] * s;
        //  multiply 'line_to_add' with 's' and subtract it from the current line 'i'
        for (AlongRowIterator iter = values_to_add.begin(); iter != values_to_add.end(); ++iter)
        {
          const double      value       = iter->second;
          const std::size_t current_col = iter->first;

          if (!value) continue;

          system_matrix(line1, current_col) -= value * s ;
        }

        system_matrix.row(line1).erase(column_to_zero);

      } // combine_lines
    } // namespace detail

    /** @brief Eliminates all odd spherical harmonics expansion coefficients in the off diagonals of S^oo from the system matrix by line operations.
     *
     * @tparam FullMatrixType   A the moment this should be the viennashe built-in matrix type (see linalg_util.hpp)
     * @tparam VectorType       Should be ublas and std::vector compatible
     * @param system_matrix     The full system matrix
     * @param rhs               The full right hand side
     * @param num_even          The number of even unkowns (needed to find S^oo)
     */
    template <typename FullMatrixType, typename VectorType>
    void diagonalise_odd2odd_coupling_matrix(FullMatrixType & system_matrix, VectorType & rhs, std::size_t num_even)
    {
      typedef typename FullMatrixType::row_type    RowType;
      typedef typename RowType::iterator           AlongRowIterator;

      const std::size_t N = rhs.size();

      // avoid round off errors by dividing through large numbers!
      //viennashe::math::row_normalize_system(system_matrix, rhs);

      // for each line in S^oo ; starting at the TOP + 1
      for ( std::size_t i = num_even + 1; i < N; ++i)
      {
        // get the line
        RowType & line = system_matrix.row(i);
        // get the off diagonal and add line
        for (AlongRowIterator cit = line.begin(); cit != line.end(); ++cit)
        {
          const std::size_t column_id = cit->first;

          if (column_id >= i)
            break; // stay left to the diagonal

          if ( column_id < num_even)
            continue; // skip elements in S^oe
          else
          {
            // column_id < i ... ok we are left to the diagonal
            // now combine lines
            const std::size_t line_to_add = i - (i - column_id);

            //std::cout << "i = " << i << " # col = " << column_id << " # dist = " << i-column_id << " # line_to_add = " << line_to_add << std::endl;

            if (line_to_add < num_even)
            {
              std::cerr << "Oh no! line_to_add < num_even " << std::endl;
              std::cerr << "Oh no! " << i << " - " << column_id << " < " << num_even << std::endl;
              throw "Oh no! line_to_add < num_even ";
            }

            viennashe::she::detail::combine_lines( system_matrix, rhs, i, line_to_add, column_id );

            break; // done :)
          }
        }// get off diagonal

      } // for each line in S^oo

      for (std::size_t i = num_even - 1; i < N; ++i)
      {
        RowType & col = system_matrix.row(i);
        for (AlongRowIterator cit = col.begin(); cit != col.end(); )
        {
         if (!cit->second)
         {
            AlongRowIterator cit_delete = cit;
            ++cit;
            col.erase(cit_delete);
          }
          else
          {
            ++cit;
          }
        } // for each column
      } // for each row

      // for each line in S^oo ; starting at the BOTTOM
      for ( std::size_t i = N - 1 /* index origin 0 */;
            i >= num_even; --i)
      {
        // get the line
        RowType & line = system_matrix.row(i);
        // get the off diagonal and add line
        for (AlongRowIterator cit = line.begin(); cit != line.end(); ++cit)
        {
          const std::size_t column_id = cit->first;

          if (column_id < i && column_id >= num_even)
          {
            std::cerr << "Oh no! There are still values left to the diagonal!" << std::endl;
            std::cerr << " column_id = " << column_id << " # i = " << i << " # num_even = " << num_even << std::endl;
            std::cerr << cit->second << std::endl;
            throw "Oh no! There are still values left to the diagonal!";
          }


          if (column_id < num_even)
          {
            continue; // skip elements in S^oe
          }
          else if (column_id == i)
          {
            continue; // skip the diagonal element
          }
          else
          {
            // column_id > i ... ok we are *right* to the diagonal
            // now combine lines
            const std::size_t line_to_add = i + (column_id - i);

            if (line_to_add >= N)
            {
              std::cerr << "Oh no! line_to_add >= N # " << line_to_add << " >= " << N << std::endl;
              throw "Oh no! line_to_add >= N ";
            }

            viennashe::she::detail::combine_lines( system_matrix, rhs, i, line_to_add, column_id );

            break; // done :)
          }
        }// get off diagonal

      } // for each line in S^oo


      // FINALLY: Heal the sparse matrix; find true zeros and eliminate
      for (std::size_t i = num_even - 1; i < N; ++i)
      {
        RowType & col = system_matrix.row(i);
        for (AlongRowIterator cit = col.begin(); cit != col.end(); )
        {
          if (!cit->second)
          {
            AlongRowIterator cit_delete = cit;
            ++cit;
            col.erase(cit_delete);
          }
          else
          {
            ++cit;
          }
        } // for each column
      } // for each row

    } // diagonalise_odd2odd_coupling_matrix



    /** @brief Eliminates all odd spherical harmonics expansion coefficients from the system matrix
     *
     * @param system_matrix     The full system matrix
     * @param rhs               The full right hand side
     * @param compressed_matrix The new matrix with only the even harmonics unknowns
     * @param compressed_rhs    The new right hand side vector with only the even harmonics unknowns
     */
    template <typename FullMatrixType, typename CompressedMatrixType, typename VectorType>
    void eliminate_odd_unknowns(FullMatrixType & system_matrix, VectorType const & rhs,
        CompressedMatrixType & compressed_matrix, VectorType & compressed_rhs)
    {
      typedef typename FullMatrixType::row_type     RowType;
      typedef typename FullMatrixType::iterator2    AlongRowIterator;

      std::size_t num_even = compressed_matrix.size1();

      for (std::size_t i=0; i<num_even; ++i)
      {
        //write new rhs entry:
        compressed_rhs[i] = rhs[i];

        RowType & row = system_matrix.row(i);
        for ( AlongRowIterator iter  = row.begin();
                               iter != row.end();
                             ++iter )
        {
          if ( iter->first < num_even ) //even unknown
          {
            compressed_matrix(i, iter->first) += iter->second;
          }
          else //odd unknown
          {
            //check for valid diagonal entry:
            double odd_diagonal_entry = -1.0;

            RowType & odd_row = system_matrix.row(iter->first);
            for ( AlongRowIterator odd_iter  = odd_row.begin();
                                   odd_iter != odd_row.end();
                                 ++odd_iter )
            {
              if (iter->first == odd_iter->first)
              {
                odd_diagonal_entry = odd_iter->second;
                break;
              }
            }

            if ( odd_diagonal_entry <= 0.0 )
            {
              log::error() << "ERROR in eliminate_odd_unknowns(): Diagonal entry " << iter->first
                  << " not positive or not existing: " << odd_diagonal_entry << std::endl;
              throw viennashe::she::invalid_matrixelement_exception("eliminate_odd_unknowns(): Diagonal entry not positive or not existing!", odd_diagonal_entry);
            }

            for ( AlongRowIterator odd_iter  = odd_row.begin();
                                   odd_iter != odd_row.end();
                                 ++odd_iter )
            {
              if (odd_iter->first < num_even)
              {
                compressed_matrix(i, odd_iter->first) -= iter->second * odd_iter->second / odd_diagonal_entry;
              }
              else if (iter->first != odd_iter->first)
              {
                log::error() << "WARNING in eliminate_odd_unknowns(): Odd2Odd block matrix has offdiagonal entry at ("
                    << iter->first << "," << odd_iter->first << ")" << std::endl;
              }
            }

            //write RHS entry (pay attention to sign flips!)
            compressed_rhs[i] -= iter->second * rhs[iter->first] / odd_diagonal_entry;
          }
        }
      }

    } //eliminate_odd_unknowns

    /** @brief Recovers the odd-order unknowns from the even-order unknowns.
     *
     * For a system matrix
     *
     * S = (S^ee, S^eo )
     *     (S^oe, S^oo ),
     *
     * the odd-order unknowns f^o are recovered from the even unknowns f^e by f^o = (S^oo)^-1 * S^oe * f^e. Note that S^oo is diagonal.
     */
    template <typename MatrixT, typename VectorT>
    VectorT recover_odd_unknowns(MatrixT const & full_matrix,
                                 VectorT const & full_rhs,
                                 VectorT const & compressed_result)
    {
      typedef typename MatrixT::row_type           RowType;
      typedef typename MatrixT::const_iterator2    AlongRowIterator;

      VectorT full_result(full_rhs.size());
      size_t num_even = compressed_result.size();

      for (std::size_t i=0; i<full_matrix.size1(); ++i)
      {
        RowType const & row_i = full_matrix.row(i);

        //write even coefficients:
        if ( i < num_even )
          full_result[i] = compressed_result[i];
        else
        {
          //compute odd coefficients:
          double diagonal_entry = 0.0;
          full_result[i] = full_rhs[i];
          for ( AlongRowIterator iter  = row_i.begin();
                                 iter != row_i.end();
                               ++iter )
          {
            if ( iter->first < num_even )
            {
              full_result[i] -= iter->second * full_result[iter->first];
            }
            else if ( i == iter->first )
            {
              diagonal_entry = iter->second;
            }
            else if( iter->second != 0.0) // check for *exactly* 0.0
            {
              log::warn() << "Non-diagonal entry found in odd-to-odd block!"  << std::endl;
              log::info<log_recover_odd_unknowns > () << "col_iter.index1() = " << i << std::endl;
              log::info<log_recover_odd_unknowns > () << "col_iter.index2() = " << iter->first << std::endl;
              log::info<log_recover_odd_unknowns > () << "*col_iter = " << iter->second << std::endl;
              assert(bool("Non-diagonal entry found in odd-to-odd block!"));
            }
          }
          assert(diagonal_entry > 0.0 && bool("FATAL ERROR: Diagonal entry missing when recovering odd coefficients!"));
          full_result[i] /= diagonal_entry;
        }
      }

      return full_result;
    } //recover_odd_unknowns


  } //namespace she
} //namespace viennashe

#endif
