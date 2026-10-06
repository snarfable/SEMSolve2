#include <basis.h>

void basis_init( basis_t *basis_in, long int dim_in, long int deg_in, long int qdeg_in )
{
	basis_in->dim = dim_in;
	basis_in->deg = deg_in;
	basis_in->pdim = binomial( dim_in + deg_in, dim_in );
	basis_in->qdeg = qdeg_in;
	basis_in->qpdim = binomial( dim_in + qdeg_in, dim_in );
}

void basis_init_std( basis_t *basis_in, long int dim_in, long int deg_in, long int qdeg_in )
{
	basis_init( basis_in, dim_in, deg_in, qdeg_in );
	basis_load_std_nodes( basis_in );
	basis_gen_coeffs( basis_in );
	basis_load_std_quad_nodes( basis_in );
	basis_gen_quad_coeffs( basis_in );
}

void basis_free( basis_t *basis_in )
{
	free( basis_in->nodes );
	free( basis_in->coeffs );
	free( basis_in->qnodes );
	free( basis_in->qcoeffs );
}

void basis_load_nodes( basis_t *basis_in, char *file_in )
{
	
}

void basis_load_std_nodes( basis_t *basis_in )
{
	long int i,j;
	long int dim = basis_in->dim;
	long int pdim = basis_in->pdim;
	long int deg = basis_in->deg;
	long int cmb[dim];
	double step = 1.0 / (double) deg;

	basis_in->nodes = (double*) malloc( dim * pdim * sizeof(double) );

	/* Make dim-combinations of dim + deg objects */
	for(i=0;i<pdim;i++)
	{
		basis_term_exponent( dim, deg, i, cmb );
		for(j=0;j<dim;j++)
			basis_in->nodes[i*dim+j] = (double) cmb[j] * step;
	}
}

void basis_load_quad_nodes( basis_t *basis_in, char *file_in )
{
	
}

void basis_load_std_quad_nodes( basis_t *basis_in )
{
	long int i,j;
	long int dim = basis_in->dim;
	long int qpdim = basis_in->qpdim;
	long int qdeg = basis_in->qdeg;
	long int cmb[dim];
	double step = 1.0 / (double) qdeg;

	basis_in->qnodes = (double*) malloc( dim * qpdim * sizeof(double) );

	/* Make dim-combinations of dim + qdeg objects */
	for(i=0;i<qpdim;i++)
	{
		basis_term_exponent( dim, qdeg, i, cmb );
		for(j=0;j<dim;j++)
			basis_in->qnodes[i*dim+j] = (double) cmb[j] * step;
	}
}

void basis_print_interp_summary( basis_t *basis_in )
{
	long int i,j,k,dim = basis_in->dim,deg = basis_in->deg,pdim = basis_in->pdim,exp[dim];
	double *coeffs = basis_in->coeffs,*nodes = basis_in->nodes;

	printf( "Interpolation Nodes:\n" );

	/* Iterate through interpolation nodes */
	for(i=0;i<pdim;i++)
	{
		printf( "node %d = ", i );
		for(j=0;j<dim;j++)
			printf( "%f ", nodes[i*dim+j] );
		printf( "\n" );
		printf( "\tinterp coeffs:\n" );
		for(j=0;j<pdim;j++)
		{
#ifdef SEMSOLVE2_BASIS_USE_LAPACK /* Meaning that the coeffs matrix is stored transposed */
			printf( "\t\t%f on ", coeffs[i*pdim+j] );
#else
			printf( "\t\t%f on ", coeffs[j*pdim+i] );
#endif
			printf( "%d = ", j );
			basis_term_exponent( dim, deg, j, exp );
			for(k=0;k<dim;k++)
				printf( "%d ", exp[k] );
			printf( "\n" );
		}
	}
}

void basis_print_quad_summary( basis_t *basis_in )
{
	long int i,j,dim = basis_in->dim,qpdim = basis_in->qpdim;
	double *qcoeffs = basis_in->qcoeffs,*qnodes = basis_in->qnodes;

	printf( "Quadrature Nodes:\n" );

	/* Iterate through quadrature nodes */
	for(i=0;i<qpdim;i++)
	{
		printf( "qnode %d = ", i );
		for(j=0;j<dim;j++)
			printf( "%f ", qnodes[i*dim+j] );
		printf( ":: coeff = %f\n", qcoeffs[i] );
	}
}

/**
 * \param dim_in The spatial dimension; generally get this from the basis_t struct
 * \param deg_in The maximum polynomial degree desired
 * \param vec_in The exponent whose index is desired
 * \return The index
 */
long int basis_term_index( long int dim_in, long int deg_in, long int *vec_in )
{
	long int i,j,dcl = 0,idx = 0;
	for(i=0;i<dim_in;i++)
		dcl += vec_in[i];
	if(dcl>deg_in)
		return -1;
	for(i=dim_in-1;i>=0;i--) // start at most significant digit
	{
		dcl = deg_in; // calculate the degree ceiling
		for(j=i+1;j<dim_in;j++)
			dcl -= vec_in[j];
		if(dcl<0)
			return -1;
		for(j=0;j<vec_in[i];j++)
			idx += basis_poly_dim( i, dcl - j );
	}
	return idx;
}

int basis_term_exponent( long int dim_in, long int deg_in, long int idx_in, long int *exp_out )
{
	long int i,j,sum,tmp,tdx = idx_in,tdg = deg_in;
	if( idx_in >= basis_poly_dim( dim_in, deg_in ) )
		return -1;
	for(i=dim_in-1;i>=0;i--)
	{
		sum = 0;
		for(j=0;;j++)
		{
			tmp = basis_poly_dim( i, tdg - j );
			if( sum + tmp > tdx )
			{
				// set exp_out[i] = j + 1;
				exp_out[i] = j;
				tdg = tdg - j;
				tdx = tdx - sum; // sum BEFORE adding tmp!!!
				break;
			}
			else
				sum += tmp;
		}
	}
	return 0;

}

/**
 * I love this function. This function does dimension-generalized
 * (i.e. for any dimension) integration of arbitrary polynomial terms.
 * The dim_in term is the dimension, and the exp_in variable is a pointer
 * to dim_in long integers defining the exponents for each independent
 * variable.
 */
double basis_term_integral_simplex( long int dim_in, long int *exp_in )
{
	long int num = 1,dnm = 0;
	long int i,j;
	for(i=0;i<dim_in;i++)
	{
		num *= factorial( exp_in[i] );
		dnm += exp_in[i];
	}
	dnm += dim_in;
	return (double) num / (double) factorial( dnm );
}

/**
 * This function assumes that the interpolation nodes have already
 * been loaded into the basis_t structure.
 * @param basis_in The basis into which to load the coefficients
 */
int basis_gen_coeffs( basis_t *basis_in )
{
	long int i,j,k,dim = basis_in->dim,pdim = basis_in->pdim,deg = basis_in->deg,exp[dim];
	double *nodes = basis_in->nodes,mat[pdim*pdim],*tmp;

	/* Allocate basis_in->coeffs */
	tmp = (double*) malloc( pdim * pdim * sizeof(double) );
	if( tmp == NULL )
		return -1;
	basis_in->coeffs = tmp;

#ifdef SEMSOLVE2_BASIS_USE_LAPACK
	int n = (int) pdim,info,ipiv[pdim];
#else
	dmatrix_t dmat,rmat;
#endif

	/* Make the righthand side unity matrix */
	for(i=0;i<pdim*pdim;i++)
	{
		if( i % pdim == i / pdim )
			basis_in->coeffs[i] = 1.0;
		else
			basis_in->coeffs[i] = 0.0;
	}

	/* Make the Vandermonde matrix */
#ifndef SEMSOLVE2_BASIS_USE_LAPACK
	for(j=0;j<pdim;j++)
	{
		/* Form the function exponents */
		basis_term_exponent( dim, deg, j, exp );
		for(i=0;i<pdim;i++)
		{
			/* Now evaluate it at interpolation node i */
			mat[i*pdim+j] = 1.0;
			for(k=0;k<dim;k++)
				mat[i*pdim+j] *= powf( nodes[i*dim+k], (double) exp[k] );
		}
	}

	/* Set up lefthand side data structure */
	dmat.start = mat;
	dmat.rsize = pdim;
	dmat.nrows = pdim;
	dmat.ncols = pdim;

	/* Set up righthand side data structure */
	rmat.start = basis_in->coeffs;
	rmat.rsize = pdim;
	rmat.nrows = pdim;
	rmat.ncols = pdim;

	/* Here, use the basic LU linear solver */
	dmatrix_lu_dcmp( &dmat );
	dmatrix_fb_sub_mrhs( &dmat, &rmat ); /* This function still needs to be written */
#else
	for(j=0;j<pdim;j++)
	{
		/* Form the function exponents */
		basis_term_exponent( dim, deg, j, exp );
		for(i=0;i<pdim;i++)
		{
			/* Now evaluate it at interpolation node i */
			mat[j*pdim+i] = 1.0;
			for(k=0;k<dim;k++)
				mat[j*pdim+i] *= powf( nodes[i*dim+k], (double) exp[k] );
		}
	}

	/* Here using lapack */
	dgesv_( &n, &n, mat, &n, ipiv, basis_in->coeffs, &n, &info );
#endif
	return 0;
}

int basis_gen_quad_coeffs( basis_t *basis_in )
{
	long int i,j,k,dim = basis_in->dim,qpdim = basis_in->qpdim,qdeg = basis_in->qdeg,exp[dim];
	double *qnodes = basis_in->qnodes,mat[qpdim*qpdim],*tmp;

	/* Allocate basis_in->coeffs */
	tmp = (double*) malloc( qpdim * sizeof(double) );
	if( tmp == NULL )
		return -1;
	basis_in->qcoeffs = tmp;

#ifdef SEMSOLVE2_BASIS_USE_LAPACK
	int n = (int) qpdim,nrhs = 1,ipiv[qpdim],info;
#else
	dmatrix_t dmat;
	dvector_t rvec;
#endif

	/* Make the righthand side unity matrix */
	for(i=0;i<qpdim;i++)
	{
		basis_term_exponent( dim, qdeg, i, exp );
		basis_in->qcoeffs[i] = basis_term_integral_simplex( dim, exp );
	}

	/* Make the Vandermonde matrix */
#ifndef SEMSOLVE2_BASIS_USE_LAPACK
	for(j=0;j<qpdim;j++)
	{
		/* Form the function exponents */
		basis_term_exponent( dim, qdeg, j, exp );
		for(i=0;i<qpdim;i++)
		{
			/* Now evaluate it at interpolation node i */
			mat[j*qpdim+i] = 1.0;
			for(k=0;k<dim;k++)
				mat[j*qpdim+i] *= powf( qnodes[i*dim+k], (double) exp[k] );
		}
	}

	/* Set up lefthand side data structure */
	dmat.start = mat;
	dmat.rsize = qpdim;
	dmat.nrows = qpdim;
	dmat.ncols = qpdim;

	/* Set up righthand side data structure */
	rvec.start = basis_in->qcoeffs;
	rvec.rsize = 1;
	rvec.length = qpdim;

	/* Here, use the basic LU linear solver */
	dmatrix_lu_dcmp( &dmat );
	dmatrix_fb_sub_vrhs( &dmat, &rvec );
#else
	for(j=0;j<qpdim;j++)
	{
		/* Form the function exponents */
		basis_term_exponent( dim, qdeg, j, exp );
		for(i=0;i<qpdim;i++)
		{
			/* Now evaluate it at interpolation node i */
			mat[i*qpdim+j] = 1.0;
			for(k=0;k<dim;k++)
				mat[i*qpdim+j] *= powf( qnodes[i*dim+k], (double) exp[k] );
		}
	}

	/* Here using lapack */
	dgesv_( &n, &nrhs, mat, &n, ipiv, basis_in->qcoeffs, &n, &info );
#endif
	return 0;
}

/**
 * This function returns the gradient of a basis function given its index
 * in the basis structure given as input and evaluated at @e loc_in.
 * @param basis_in The basis
 * @param findex_in The function index
 * @param loc_in The point in the standard simplex at which to evaluate the gradient
 * @param grad_out The location into which to deposit the gradient value
 * @return The value is 0 if no error, otherwise something else (?)
 */
int basis_func_gradient( basis_t *basis_in, long int fidx_in, double *loc_in, double *grad_out )
{
	long int i,j,k,exp[basis_in->dim],dim = basis_in->dim,pdim = basis_in->pdim,deg = basis_in->deg;
	double prd;

	/* For each variable with which to differentiate */
	for(i=0;i<dim;i++)
	{
		/* For each polynomial term in function fidx_in */
		grad_out[i] = 0.0;
		for(j=0;j<pdim;j++)
		{
			basis_term_exponent( dim, deg, j, exp );
			if( exp[i] > 0 )
				prd = (double) exp[i] * powf( loc_in[i], (double) ( exp[i] - 1 ) ); /* This is the term to differentiate */
			else /* Since the exponent is 0, then the derivative vanishes */
				prd = 0.0;
			for(k=0;k<dim;k++)
				if( k != i )
					prd *= powf( loc_in[k], (double) exp[k] );
#ifdef SEMSOLVE2_BASIS_USE_LAPACK
			grad_out[i] += basis_in->coeffs[fidx_in*pdim+j] * prd;
#else
			grad_out[i] += basis_in->coeffs[j*pdim+fidx_in] * prd;
#endif
		}
	}

	return 0;
}

/**
 * This function just calculates the number of arrangements of deg_in indistinguishable
 * balls into dim_in distinguishable cells.
 * @param dim_in The spatial dimension
 * @param deg_in The maximum degree of the polynomial
 * @return binomial( dim_in + deg_in, dim_in )
 */
long int basis_poly_dim( long int dim_in, long int deg_in )
{
	// returns C( dim_in + deg_in, dim_in )
	return binomial( dim_in + deg_in, dim_in );
}

double basis_plane_test( long int dim_in, long int nplpts_in, double *simpx_in, long int *plane_in, double *point_in )
{
	int i,j;
	int n,nrhs;
	int ipiv[dim_in+1];
	double sum;
	double vdm[dim_in*dim_in+2*dim_in+1];
	double rhs[dim_in+1];

	// make the Vandermonde matrix
	for(i=0;i<dim_in+1;i++) // column-major -> row = function = component, column = point
		for(j=0;j<dim_in+1;j++)
			vdm[(dim_in+1)*i+j] = ( i == 0 )?( 1.0 ):( simpx_in[j*dim_in+i-1] );

	// make the righthand side
	for(i=0;i<dim_in+1;i++)
		rhs[i] = 1.0;
	for(i=0;i<nplpts_in;i++)
		rhs[plane_in[i]] = 0.0;

	n = (int) dim_in + 1;
	nrhs = 1;
	dgesv_( &n, &nrhs, vdm, &n, ipiv, rhs, &n, &i );
	sum = 0.0;
	for(i=0;i<dim_in+1;i++)
		sum += ( i == 0 )?( rhs[i] ):( rhs[i] * point_in[i-1] );
	return sum;
}

/**
 * This function converts an input point, point_in, from cartesian coordinates to
 * barycentric coordinates based on an input simplex in simp_in.
 */
void basis_barycentric( long int dim_in, double *point_in, double *simp_in, double *bary_out )
{
	long int i,j;
	double mat[dim_in*dim_in];
	int n = (int)dim_in,nrhs = 1,ipiv[dim_in],info;

	/* Making the T matrix */
	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			mat[j*dim_in+i] = simp_in[j*dim_in+i] - simp_in[dim_in*dim_in+i];
	for(i=0;i<dim_in;i++)
		bary_out[i] = point_in[i] - simp_in[dim_in*dim_in+i];
	dgesv_( &n, &nrhs, mat, &n, ipiv, bary_out, &n, &info );
	bary_out[dim_in] = 1.0;
	for(i=0;i<dim_in;i++)
		bary_out[dim_in] = bary_out[dim_in] - bary_out[i];
}

double basis_eval_interp_func( basis_t *basis_in, long int fidx_in, double *point_in )
{
	long int i,j,dim = basis_in->dim,pdim = basis_in->pdim,exp[dim],deg = basis_in->deg;
	double prd,val;

	val = 0.0;
	for(i=0;i<pdim;i++)
	{
		basis_term_exponent( dim, deg, i, exp );
		prd = 1.0;
		for(j=0;j<dim;j++)
			prd *= powf( point_in[j], (double) exp[j] );
#ifdef SEMSOLVE2_BASIS_USE_LAPACK
		val += basis_in->coeffs[fidx_in*pdim+i] * prd;
#else
		val += basis_in->coeffs[i*pdim+fidx_in] * prd;
#endif
	}
	return val;
}
