/* ============================================================================
 Copyright (c) 2011-2019, Institute for Microelectronics,
 Institute for Analysis and Scientific Computing,
 TU Wien.

 -----------------
 ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
 -----------------

 http://viennashe.sourceforge.net/

 Lead Developer:  Karl Rupp                           rupp@iue.tuwien.ac.at
 Developer:  Felipe Senra Riberio                ribeiro@iue.tuwien.ac.at

 License:         MIT (X11), see file LICENSE in the base directory
 =============================================================================== */

#ifndef PETSC_HPP_
#define PETSC_HPP_

#include <petscksp.h>
#include <mpi.h>
#include <iterator>
#include <vector>
#include <numeric>
#include <memory>
#include <unordered_map>

namespace viennashe
{
  namespace solvers
  {
    /** @brief Class Encapsulate The PETSc COMM WORLD
     *        -Should be singleton.
     *  **/
    template<typename NumericT, typename VectorType>
      class PETSC_vector;
    template<typename NumericT, typename VectorType>
      class PETSC_matrix;

    template<typename NumericT, typename VectorType>
      class PETSC_handler
      {
      private:

	PETSC_handler (int argc, char **argv) :
	    comm (MPI_COMM_WORLD), self (PETSC_COMM_SELF)
	{
	  MPI_Comm_size (PETSC_COMM_WORLD, &sizeW);
	  MPI_Comm_rank (PETSC_COMM_WORLD, &rankW);
	  MPI_Comm_size (PETSC_COMM_SELF, &size);
	  MPI_Comm_rank (PETSC_COMM_SELF, &rank);
	}
	void
	operator= (PETSC_handler const&); // Don't implement
	static std::unique_ptr<PETSC_handler> singleton;
	MPI_Comm comm, self;
	PetscViewer viewer;
	PetscMPIInt rankW, sizeW, rank, size;

      public:
	// Implement methods for setting the PETSC env.
	static PETSC_handler&
	get_handler ()
	{
	  if (singleton == nullptr)
	    singleton.reset (new PETSC_handler);
	  return *singleton;
	}

	static PETSC_handler&
	get_handler (int argc, char **argv)
	{
	  if (singleton == nullptr)
	    singleton.reset (new PETSC_handler (argc, argv));
	  return *singleton;
	}
	PETSC_handler (PETSC_handler const&) = delete;

	MPI_Comm
	getComm () const
	{
	  return comm;
	}

	void
	setComm (MPI_Comm comm)
	{
	  this->comm = comm;
	}

	PetscMPIInt
	getRank () const
	{
	  return rank;
	}

	PetscMPIInt
	getRankW () const
	{
	  return rankW;
	}

	MPI_Comm
	getSelf () const
	{
	  return self;
	}

	PetscMPIInt
	getSize () const
	{
	  return size;
	}

	PetscMPIInt
	getSizeW () const
	{
	  return sizeW;
	}
      }
      ;

    template<typename NumericT, typename VectorType>
      std::unique_ptr<PETSC_handler<NumericT, VectorType>> PETSC_handler<
	  NumericT, VectorType>::singleton = nullptr;

    /** @brief Wrapper class to the PETSC matrix
     *  **/
    template<typename NumericT, typename VectorType>
      class PETSC_matrix
      {
	typedef PETSC_matrix<NumericT, VectorType> self_type;
      public:

	typedef std::size_t size_type;
	typedef size_type SizeType;

      public:

	/** @brief Wrapper class to the PETSC matrix
	 *  @param copymat - Matrix to be copied to PETSc
	 *  **/
	PETSC_matrix (const viennashe::math::sparse_matrix<NumericT> &copymat) :
	      comm (PETSC_COMM_WORLD), self (
		PETSC_COMM_SELF),its (0),i (0), j (0), index (0), Creason (KSP_CONVERGED_RTOL_NORMAL), iter (0)
	{
	  typedef typename viennashe::math::sparse_matrix<NumericT>::RowType RowType;
	  // typedef typename RowType::iterator iterator2;
	  typedef typename RowType::const_iterator const_iterator2;

	  //* Get some information from the matrix*/
	  MPI_Comm_size (PETSC_COMM_WORLD, &size);
	  MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
	  int n = copymat.size1 ();
	  int m = copymat.size2 ();
	  long unsigned int local = copymat.size2();
	  rstart = 0, rend =m ;
	  size1_ = n;
	  size2_ = m;
	  long unsigned int nnz = 0;
	  for (unsigned long int k = 0; k < copymat.size1(); k++)
  	    {
  	     nnz = nnz < copymat.row(k).size()? copymat.row(k).size() :nnz;
  	    }
	  nnz++;
	  PetscInt *col = new PetscInt[nnz];
  	  nghosts = 0;
  	  double *value = new double[nnz];
  	  this->index = 0;
	  ierr = MatCreate (PETSC_COMM_WORLD, &Apetsc);
	  CHKERRV (ierr);

	  //* Split matrix for SHE only//
	  if(PETSC_matrix::ninstan == 2){
	    long unsigned int initial;
	    long unsigned int total;
	    RowType const &Row = copymat.row (local-1);
	    offset = 0;
	    if(rank != size-1) offset = Row.rbegin()->first -Row.begin()->first ;
	    local = copymat.size2();
	    MPI_Scan(&local, &initial, 1, MPI_UNSIGNED_LONG , MPI_SUM, PETSC_COMM_WORLD);
	    MPI_Allreduce(&local, &total, 1, MPI_UNSIGNED_LONG , MPI_SUM, PETSC_COMM_WORLD);
	    PETSC_matrix::ninstan = -1;
	    PetscBarrier (NULL);
	    m = n=size2_ = total;
	    rstart = initial-local;
	    rend = initial;
	    ierr = MatSetSizes (Apetsc,local, local,  PETSC_DECIDE, PETSC_DECIDE);
	    CHKERRV (ierr);

	  }else{
	      ierr = MatSetSizes (Apetsc, PETSC_DECIDE, PETSC_DECIDE,  m, n);
	      offset = 0;
	      CHKERRV (ierr);
	  }

	  ierr = MatSetFromOptions (Apetsc);
	  CHKERRV (ierr);
	  ierr = MatSeqAIJSetPreallocation (Apetsc, nnz, PETSC_NULL);
	  	  CHKERRV (ierr);
	  ierr = MatMPIAIJSetPreallocation (Apetsc, nnz, PETSC_NULL, nnz,
					    PETSC_NULL);
	  CHKERRV (ierr);

	  ierr = MatSetUp (Apetsc);
	  CHKERRV (ierr);
	  ierr = MatGetLocalSize (Apetsc, &nlocal1, &nlocal2);
	  CHKERRV (ierr);
	  ierr = MatSetOption (Apetsc, MAT_NEW_NONZERO_ALLOCATION_ERR,
			       PETSC_FALSE);
	  PetscInt  nlocal;
	  for (PetscInt k = 0; k < (PetscInt)(copymat.size1()); k++)
	    {
	      RowType const &Row = copymat.row (k);
	      long nnzt = 0;

	      for (const_iterator2 it = Row.begin (); it != Row.end (); it++)
		{
		  col[nnzt] = it->first+rstart;
		  value[nnzt] = it->second;
		  nlocal = k+rstart;
		  nnzt++;
		}
	      ierr = MatSetValues (Apetsc, 1, &(nlocal), nnzt, col, value, INSERT_VALUES);
	      CHKERRV (ierr);
	    }

	  delete[] col;
	  delete[] value;
	  ierr = MatAssemblyBegin (Apetsc, MAT_FINAL_ASSEMBLY);
	  CHKERRV (ierr);
	  ierr = MatAssemblyEnd (Apetsc, MAT_FINAL_ASSEMBLY);
	  CHKERRV (ierr);
	  nconsiter++;
	  KSPCreate (comm, &ksp);
	  KSPSetOperators (ksp, Apetsc, Apetsc);
	  PETSC_matrix::ninstan++;
	}

	~PETSC_matrix ()
	{
	  MatDestroy (&Apetsc);
	  KSPDestroy (&ksp);
	}

	/** @brief Wrapper class to the PETSC matrix
	 *  @param b - right hand term
	 *  **/
	void
	solve (PETSC_vector<NumericT, VectorType> &b)
	{
	  PetscErrorCode ierr;
	  KSPSetFromOptions (ksp);
	  ////printf ("\n Solving");

	  ierr = KSPSolve (ksp, b.getVec (), b.getSolution ());CHKERRV (ierr);
	  iter++;
	  KSPGetIterationNumber (ksp, &its);
	  b.GenerateGlobal ();
	}

	size_type
	size1 () const
	{
	  return size1_;
	}
	size_type
	size2 () const
	{
	  return size2_;
	}



	typedef PETSC_vector<NumericT, VectorType> row_type;
	typedef PETSC_vector<NumericT, VectorType> RowType;

	self_type&
	operator = (double other) // note: argument passed by value!
	{
	  ierr = MatSetValue (Apetsc, i, j, other, INSERT_VALUES);
	  return *this;
	}

	PetscReal
	getNorm () const
	{
	  return norm;
	}
	KSPConvergedReason
	get_Creason ()
	{
	  KSPGetConvergedReason (ksp, &Creason);
	  return Creason;
	}

	int
	get_Citer ()
	{
	  KSPGetIterationNumber (ksp, &its);
	  return its;
	}
	void
	setNorm (PetscReal norm)
	{
	  this->norm = norm;
	}

	size_type
	getSize1 () const
	{
	  return size1_;
	}

	void
	setSize1 (size_type size1)
	{
	  size1_ = size1;
	}

	size_type
	getSize2 () const
	{
	  return size2_;
	}

	void
	setSize2 (size_type size2)
	{
	  size2_ = size2;
	}

	const Mat&
	getA ()
	{
	  return Apetsc;
	}

	int
	getIndex () const
	{
	  return index;
	}

	void
	setIndex (int index)
	{
	  this->index = index;
	}

	PetscInt
	getNlocal1 () const
	{
	  return nlocal1;
	}

	void
	setNlocal1 (PetscInt nlocal1)
	{
	  this->nlocal1 = nlocal1;
	}

	PetscInt
	getNlocal2 () const
	{
	  return nlocal2;
	}

	void
	setNlocal2 (PetscInt nlocal2)
	{
	  this->nlocal2 = nlocal2;
	}

	int
	getIter () const
	{
	  return iter;
	}

	void
	setIter (int iter)
	{
	  this->iter = iter;
	}

	PetscInt
	getRend () const
	{
	  return rend;
	}

	void
	setRend (PetscInt rend)
	{
	  this->rend = rend;
	}

	PetscInt
	getRstart () const
	{
	  return rstart;
	}

	void
	setRstart (PetscInt rstart)
	{
	  this->rstart = rstart;
	}

      private:
	PetscMPIInt rank, size;
	PetscInt rstart , rend, offset ;
	MPI_Comm comm, self;
	IS is;
	MatPartitioning part;
	PetscInt its, nlocal1, nlocal2;
	long i, j, index;
	Mat Apetsc;
	KSP ksp;
	// PC pc;
	PetscReal norm;
	KSPConvergedReason Creason;
	PetscErrorCode ierr = 0;
	size_type size1_;
	size_type size2_;
	int nghosts, iter;
	static int ninstan;
	static int nconsiter;

      };
    /** @brief Wrapper class to the PETSC vector
     *  **/
    template<typename NumericT, typename VectorType>
      class PETSC_vector
      {
	typedef PETSC_vector<NumericT, VectorType> self_type;
      public:
	typedef std::size_t size_type;
	typedef NumericT value_type;

	/** @brief Constructor to the PETSC vector Right Hand and Solution
	 *      M  Matrix of the Linear System
	 *  **/
	PETSC_vector (const std::vector<NumericT> &V,
		      PETSC_matrix<NumericT, VectorType> &M) :
	    i(0),nlocal(0), residualb (false), globalb (false),size_n (V.size ()), dad (&M)
	{
	  PetscInt  nlocal;
	  ierr = MatCreateVecs ((M.getA ()), &vec, &solution);
	  VecDuplicate (vec, &residual);
	  ierr = PetscObjectSetName ((PetscObject) vec, "RHS");
	  ierr = PetscObjectSetName ((PetscObject) solution, "Solution");
	  VecGetOwnershipRange (vec, &rstart, &rend);
	  VecGetLocalSize (vec, &nlocal);
	  sizeg = M.size2();
	  rstart = M.getRstart();
	  rend = M.getRend();
	  for (size_t k = 0; k < V.size(); k++)
	    {
	      if (V[k] != 0.0){
		nnz++;
	      }
	      nlocal = k+rstart;
	      VecSetValues (vec, 1, &(nlocal), &V[k], INSERT_VALUES);
	    }
	  VecAssemblyBegin (vec);
	  VecAssemblyEnd (vec);
	  VecAssemblyBegin (solution);
	  VecAssemblyEnd (solution);
	  VecAssemblyBegin (residual);
	  VecAssemblyEnd (residual);
	}
	~PETSC_vector ()
	{
	  VecDestroy (&vec);
	  VecDestroy (&solution);
	  VecDestroy (&residual);
	}

	/** @brief Cast from the PETSC vector to std::vector
	 *
	 *  **/
	operator std::vector<NumericT> ()    //      const
	{
	  std::vector<NumericT> result (size_n); // static or not
	  PetscScalar *array;
	  PetscInt size;
	  if (!globalb)
	    {
	      VecGetLocalSize (solution, &size);
	      if (residualb == false)
		{
		  VecGetArray (solution, &array);
		  std::vector<NumericT> result (array, array + (size));
		  VecRestoreArray (solution, &array);
		  VecDestroy (&residual);
		  VecDestroy (&solution);
		  VecDestroy (&vec);
		  return result;
		}
	      else
		{
		  VecGetLocalSize (residual, &size);
		  MatMult (dad->getA (), vec, residual);   // y <- Ax
		  VecAXPY (residual, -1., vec);  // r <- r - b
		  VecGetArray (residual, &array);
		  std::vector<NumericT> result (array, array + (size));
		  VecRestoreArray (residual, &array);
		  return result;
		}
	    }
	  else
	    {
	      VecGetLocalSize (global, &size);
	      if (size != 0)
		{
		  if ((residualb == false) && (size != 0))
		    {
		      VecGetArray (global, &array);
		      std::vector<NumericT> result (array+rstart,array+rend);
		      VecRestoreArray (global, &array);
		      VecScatterDestroy (&scat);
		      VecDestroy (&global);
		      VecDestroy (&residual);
		      VecDestroy (&solution);
		      VecDestroy (&vec);
		      return result;
		    }
		  else
		    {
		      VecGetLocalSize (residual, &size);
		      MatMult (dad->getA (), vec, residual); // y <- Ax
		      VecAXPY (residual, -1., vec);  // r <- r - b
		      VecGetArray (residual, &array);
		      std::vector<NumericT> result (array, array + (size));
		      VecRestoreArray (residual, &array);
		      return result;
		    }
		}

	    }
	  return result;
	}
//        /** @brief Send the solution to the first process( process 0)
//         *  **/
	void
	GenerateGlobal ()
	{
	  globalb = true;
	  VecScatterCreateToAll (solution, &scat, &global);
	  ierr = PetscObjectSetName ((PetscObject) global, "Global");
	  VecScatterBegin (scat, solution, global, INSERT_VALUES,
			   SCATTER_FORWARD);
	  VecScatterEnd (scat, solution, global, INSERT_VALUES,
			 SCATTER_FORWARD);
	}

	Vec
	getSolution () const
	{
	  return solution;
	}

	void
	setSolution (Vec solution)
	{
	  this->solution = solution;
	}
	operator Vec () const
	{
	  return vec;
	}

	size_type
	size () const
	{
	  return this->size_n;
	}

	size_type
	sizeLocal () const
	{
	  return this->nlocal;
	}

	Vec
	getVec () const
	{
	  return vec;
	}

	void
	setVec (Vec vec)
	{
	  this->vec = vec;
	}

      private:
	int i;
	int nnz, nlocal;
	bool residualb, globalb;
	PetscInt rstart , rend ;
	IS indexG, indexL;
	Vec vec, solution, residual, global;
	VecScatter scatter;
	size_type size_n, sizeg;
	VecScatter scat;
	ISLocalToGlobalMapping ltog;
	PETSC_matrix<NumericT, VectorType> *dad;
	PetscErrorCode ierr;

      };

    template<typename NumericT, typename VectorType>
      int PETSC_matrix<NumericT, VectorType>::ninstan = 0;
    template<typename NumericT, typename VectorType>
      int PETSC_matrix<NumericT, VectorType>::nconsiter = 0;
  }
}
#endif /* PETSC_HPP_ */
