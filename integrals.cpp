#include <integrals.h>
#include <potential.h>

/**
 * This function calculates the inner product of two interpolation functions on
 * the simplex by using the internal quadrature rule generated when the basis_init
 * was called on @e basis_in.
 * @param basis_in The basis to use
 * @param fidx1_in The index of the first function
 * @param fidx2_in The index of the second function
 * @return Returns the value of the inner product
 */
double inner_prod( basis_t *basis_in, long int fidx1_in, long int fidx2_in )
{
	long int i,j,dim = basis_in->dim,qpdim = basis_in->qpdim;
	double tot,db1,db2,*qnodes = basis_in->qnodes,*qcoeffs = basis_in->qcoeffs;

	tot = 0.0;
	for(i=0;i<qpdim;i++)
	{
		db1 = basis_eval_interp_func( basis_in, fidx1_in, qnodes + i * dim );
		db2 = basis_eval_interp_func( basis_in, fidx2_in, qnodes + i * dim );
		tot += qcoeffs[i] * db1 * db2;
	}
	return tot;
}

inline void matrix_multiply_slow( long int m, long int n, long int k, double *mat1_in, double *mat2_in, double *mat_out )
{
	long int i,j,p;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			mat_out[i*n+j] = 0.0;
			for(p=0;p<k;p++)
				mat_out[i*n+j] += mat1_in[i*k+p] * mat2_in[p*n+j];
		}
	}
}

inline double dot_product( long int n, double *vec1_in, double *vec2_in )
{
	long int i;
	double tot = 0.0;

	for(i=0;i<n;i++)
		tot += vec1_in[i] * vec2_in[i];

	return tot;
}

/**
 * IMPORTANT: This function returns an integral value which must be scaled
 * by the volume of the simplex. This can be done by multiplying by the determinant
 * of the linear coordinate transformation (the affine transformation) and then
 * dividing by the factorial of the dimension of the simplex. Notice that the determinant
 * is constant with respect to integration. That is, it does not change with change in
 * the input point during integration. Because of this, it need not be calculated
 * at each quadrature node.
 * @param basis_in The basis
 * @param fidx1_in The index of the first interpolation function
 * @param fidx2_in The index of the second interpolation function
 * @param simpx_in The global Euclidean coordinates of the simplex (i.e. n + 1 nodes for @e dim = @e n)
 */
double kinetic_inner_prod( basis_t *basis_in, long int fidx1_in, long int fidx2_in, double *simpx_in )
{
	long int i,j,dim = basis_in->dim,qpdim = basis_in->qpdim;
	double tot,db1[dim],db2[dim],op1[dim],op2[dim],*qnodes = basis_in->qnodes,*qcoeffs = basis_in->qcoeffs;
	double trn[dim*dim],inv[dim*dim];

#ifdef SEMSOLVE2_BASIS_USE_LAPACK
	int n,info,ipiv[dim],rhs = 1;
	double alpha = 1.0,beta = 0.0;
	char t = 'n';
#else
	dmatrix_t dmat,imat;
#endif

	/* Make the affine transformation matrix from local to global */
	for(i=0;i<dim*dim;i++)
	{
		if( i % dim == i / dim )
			inv[i] = 1.0;
		else
			inv[i] = 0.0;
	}

#ifdef SEMSOLVE2_BASIS_USE_LAPACK
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			trn[i*dim+j] = simpx_in[(j+1)*dim+i] - simpx_in[i];

	/* Invert using lapack */
	n = (int) dim;
	dgesv_( &n, &n, trn, &n, ipiv, inv, &n, &info );
#else
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			trn[j*dim+i] = simpx_in[(j+1)*dim+i] - simpx_in[i];

	/* The affine matrix */
	dmat.start = trn;
	dmat.rsize = dim;
	dmat.nrows = dim;
	dmat.ncols = dim;

	/* The output */
	imat.start = inv;
	imat.rsize = dim;
	imat.nrows = dim;
	imat.ncols = dim;

	/* Invert using my own personal garbage bag of linear algebra functions */
	dmatrix_lu_dcmp( &dmat );
	dmatrix_fb_sub_mrhs( &dmat, &imat );
#endif

	/* Set inner sum to zero */
	tot = 0.0;

	/* Evaluate both gradients at each quadrature node */
	for(i=0;i<qpdim;i++)
	{
		basis_func_gradient( basis_in, fidx1_in, qnodes + i * dim, db1 );
		basis_func_gradient( basis_in, fidx2_in, qnodes + i * dim, db2 );
#ifdef SEMSOLVE2_BASIS_USE_LAPACK
		dgemm_( &t, &t, &n, &rhs, &n, &alpha, inv, &n, db1, &n, &beta, op1, &n );
		dgemm_( &t, &t, &n, &rhs, &n, &alpha, inv, &n, db2, &n, &beta, op2, &n );
#else
		matrix_multiply_slow( dim, 1, dim, inv, db1, op1 );
		matrix_multiply_slow( dim, 1, dim, inv, db2, op2 );
#endif
		tot += qcoeffs[i] * dot_product( dim, op1, op2 );
	}

	return tot;
}
