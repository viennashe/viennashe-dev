#ifndef VIENNASHE_UTIL_BLOCK_MATRIX_WRITER_HPP
#define VIENNASHE_UTIL_BLOCK_MATRIX_WRITER_HPP

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

#include <cassert>

// viennashe
#include "viennashe/math/spherical_harmonics.hpp"
#include "viennashe/physics/constants.hpp"

/** @file viennashe/util/block_matrix_writer.hpp
    @brief Writes a possibly scaled block matrix to the system matrix. Includes similar functionality for computing residual contributions.
*/

namespace viennashe
{
  namespace util
  {

    /** @brief Writes a sub-matrix of a block matrix 'coupling_matrix' into the global system matrix, scaled by 'prefactor'
    *
    * @param system_matrix    System matrix, to which the entries are written
    * @param row_index        Starting row index in system matrix
    * @param col_index        Starting column index in system matrix
    * @param prefactor        The block matrix 'coupling_matrix' is multiplied by prefactor before the write operation is carried out
    * @param coupling_matrix  The block matrix
    * @param row_iter         Iterator over row indices of 'coupling matrix'
    * @param col_iter_init    Iterator over column indices of 'coupling matrix'
    * @param elastic_scattering_roundoff_error_prevention  Flag which specifies whether the first element in the coupling matrix should be forced to zero to avoid round-off errors
    */
    template <typename SystemMatrixType,
              typename BlockMatrixType,
              typename RowIndexIterator,
              typename ColumnIndexIterator>
    void add_block_matrix(SystemMatrixType & system_matrix,
                          std::size_t row_index, std::size_t col_index,
                          double prefactor,
                          BlockMatrixType const & coupling_matrix,
                          RowIndexIterator row_iter,
                          ColumnIndexIterator const & col_iter_init,
                          bool elastic_scattering_roundoff_error_prevention = false)
    {
      //iterate over harmonics in row and column direction

      assert(prefactor == prefactor && bool("Writing nan to matrix!"));

      std::size_t row = row_index;
      for (; row_iter.valid(); ++row_iter, ++row)
      {
        assert( *row_iter < static_cast<long>(coupling_matrix.size1() ));
        std::size_t col = col_index;

        for (ColumnIndexIterator col_iter(col_iter_init);
             col_iter.valid();
             ++col_iter, ++col)
        {
          assert( *col_iter < static_cast<long>(coupling_matrix.size2() ));
          double entry = coupling_matrix(std::size_t(*row_iter), std::size_t(*col_iter));

          if (entry) //preserve sparsity structure of coupling matrices. Double-inequality to avoid warnings regarding equality to zero
          {
            if (elastic_scattering_roundoff_error_prevention && row == row_index && col == col_index) // avoid round-off errors for scattering mechanisms
              continue;
            system_matrix(row, col) += prefactor * entry;
          }
        }
      }
    }


    /** @brief Writes a sub-matrix of a block matrix 'coupling_matrix' into the global system matrix, scaled by 'prefactor'. The columns are 'folded' (summed, i.e. a matrix-vector product) according to the vector provided.
    *
    * @param system_matrix    System matrix, to which the entries are written
    * @param row_index        Starting row index in system matrix
    * @param col_index        Starting column index in system matrix
    * @param prefactor        The block matrix 'coupling_matrix' is multiplied by prefactor before the write operation is carried out
    * @param coupling_matrix  The block matrix
    * @param row_iter         Iterator over row indices of 'coupling matrix'
    * @param col_iter_init    Iterator over column indices of 'coupling matrix'
    * @param fold_vector      The vector with which the fold-operation (matrix-vector product on the columns selected by the column-iterator) is carried out.
    * @param fold_index_start Start (offset) index within the fold_vector
    */
    template <typename SystemMatrixType,
              typename BlockMatrixType,
              typename RowIndexIterator,
              typename ColumnIndexIterator,
              typename FoldVectorType>
    void add_folded_block_matrix(SystemMatrixType & system_matrix,
                                 long row_index, long col_index,
                                 double prefactor,
                                 BlockMatrixType const & coupling_matrix,
                                 RowIndexIterator row_iter,
                                 ColumnIndexIterator const & col_iter_init,
                                 FoldVectorType const & fold_vector, long fold_index_start
                                )
    {
      //iterate over harmonics in row and column direction
      assert(row_index >= 0 && "Row index must be non-negative!");
      assert(col_index >= 0 && "Col index must be non-negative!");

      assert(prefactor == prefactor && bool("Writing nan to matrix!"));

      long row = row_index;
      for (; row_iter.valid(); ++row_iter)
      {
        assert( *row_iter < static_cast<long>(coupling_matrix.size1()) );

        double folded_val = 0;
        std::size_t inner_iters = 0;
        for (ColumnIndexIterator col_iter(col_iter_init);
             col_iter.valid();
             ++col_iter, ++inner_iters)
        {
          assert( *col_iter < static_cast<long>(coupling_matrix.size2()) );
          folded_val += coupling_matrix(*row_iter, *col_iter) * fold_vector[fold_index_start + inner_iters];
        }

        if (folded_val) //preserve sparsity structure of coupling matrices. Double-inequality to avoid warnings regarding equality to zero
          system_matrix(row, col_index) += prefactor * folded_val;

        ++row;
      }
    }


    /** @brief Writes a sub-matrix of a block matrix 'coupling_matrix' to the residual vector.  The columns are 'folded' (summed, i.e. a matrix-vector product) according to the vector provided.
    *
    * @param residual         Residual vector to which the entries are written
    * @param row_index        Starting row index in system matrix
    * @param prefactor        The block matrix 'coupling_matrix' is multiplied by prefactor before the write operation is carried out
    * @param coupling_matrix  The block matrix
    * @param row_iter         Iterator over row indices of 'coupling matrix'
    * @param col_iter_init    Iterator over column indices of 'coupling matrix'
    * @param fold_vector      The vector with which the fold-operation (matrix-vector product on the columns selected by the column-iterator) is carried out.
    */
    template <typename VectorType,
              typename BlockMatrixType,
              typename RowIndexIterator,
              typename ColumnIndexIterator>
    void subtract_folded_block_vector(VectorType & residual,
                                      std::size_t row_index,
                                      double prefactor,
                                      BlockMatrixType const & coupling_matrix,
                                      RowIndexIterator row_iter,
                                      ColumnIndexIterator const & col_iter_init,
                                      double const * fold_vector)
    {
      assert(prefactor == prefactor && bool("Writing nan to vector!"));

      //iterate over harmonics in row and column direction

      std::size_t row = row_index;
      for (; row_iter.valid(); ++row_iter)
      {
        assert( *row_iter < static_cast<long>(coupling_matrix.size1() ));

        double folded_val = 0;
        std::size_t inner_iters = 0;
        for (ColumnIndexIterator col_iter(col_iter_init);
             col_iter.valid();
             ++col_iter, ++inner_iters)
        {
          assert( *col_iter < static_cast<long>(coupling_matrix.size2() ));
          folded_val += coupling_matrix(std::size_t(*row_iter), std::size_t(*col_iter)) * fold_vector[inner_iters];
        }

        if (folded_val) //preserve sparsity structure of coupling matrices. Double-inequality to avoid warnings regarding equality to zero
          residual[row] -= prefactor * folded_val;

        ++row;
      }
    }

  } //namespace util
} //namespace viennashe

#endif
