#include <slinalg.h>

void smatrix_dspmv( smatrix_t *mat_in, dvector_t *vec_in, dvector_t *vec_out )
{
	long int i,j;

	for(i=0;i<mat_in->nrows;i++)
	{
		vec_out->start[i*vec_out->rsize] = 0.0;
		for(j=0;j<mat_in->lvec[i];j++)
			vec_out->start[i*vec_out->rsize] += vec_in->start[mat_in->jvec[i][j]*vec_in->rsize] * mat_in->mvec[i][j];
	}
}

void smatrix_dspmv_sparse( smatrix_t *mat_in, smatrix_t *vec_in, long int col1_in, smatrix_t *vec_out, long int col2_in )
{
	long int i,j;
	double sum,tmp;
	for(i=0;i<mat_in->nrows;i++)
	{
		sum = 0.0;
		for(j=0;j<mat_in->lvec[i];j++)
		{
			smatrix_entry( vec_in, mat_in->jvec[i][j], col1_in, &tmp );
			mat_in->mvec[i][j];
			sum += tmp * mat_in->mvec[i][j];
		}
		smatrix_insert( vec_out, i, col2_in, sum );
	}
}

void smatrix_dspmv_sparse_dvector( smatrix_t *mat_in, smatrix_t *vec_in, long int col1_in, dvector_t *vec_out )
{
	long int i,j;
	double sum,tmp;
	for(i=0;i<mat_in->nrows;i++)
	{
		sum = 0.0;
		for(j=0;j<mat_in->lvec[i];j++)
		{
			smatrix_entry( vec_in, mat_in->jvec[i][j], col1_in, &tmp );
			mat_in->mvec[i][j];
			sum += tmp * mat_in->mvec[i][j];
		}
		vec_out->start[i*vec_out->rsize] = sum;
	}
}

void smatrix_lin_com( double cof1_in, smatrix_t *mat1_in, long int col1_in, double cof2_in, smatrix_t *mat2_in, long int col2_in, smatrix_t *mat_out, long int col_out )
{
	long int i,j;
	double ent1,ent2,tmp;
	assert( mat1_in->nrows == mat2_in->nrows && mat2_in->nrows == mat_out->nrows );
	for(i=0;i<mat1_in->nrows;i++)
	{
		smatrix_entry( mat1_in, i, col1_in, &ent1 );
		smatrix_entry( mat2_in, i, col2_in, &ent2 );
		tmp = cof1_in * ent1 + cof2_in * ent2;
		smatrix_insert( mat_out, i, col_out, tmp );
	}
}

double smatrix_col_dot_prod( smatrix_t *mat1_in, long int col1_in, smatrix_t *mat2_in, long int col2_in )
{
	long int i,j;
	double tmp1,tmp2,sum;
	assert( mat1_in->nrows == mat2_in->nrows );
	sum = 0.0;
	for(i=0;i<mat1_in->nrows;i++)
	{
		smatrix_entry( mat1_in, i, col1_in, &tmp1 );
		smatrix_entry( mat2_in, i, col2_in, &tmp2 );
		sum += tmp1 * tmp2;
	}
	return sum;
}

double smatrix_col_dot_prod_dvector( smatrix_t *mat_in, long int col_in, dvector_t *vec_in )
{
	long int i,j;
	double tmp,sum;
	assert( mat_in->nrows == vec_in->length );
	sum = 0.0;
	for(i=0;i<mat_in->nrows;i++)
	{
		smatrix_entry( mat_in, i, col_in, &tmp );
		sum += tmp * vec_in->start[i*vec_in->rsize];
	}
	return sum;
}

/**
 * This is a sparse LU decomposition algorithm.
 * \param mat_in The input matrix to decompose into lower and upper triangular components
 * \param mat_out The lower and upper triangular decomposition products are both stored in \c mat_out
 *	with the exception of the diagonal of the upper triangular component, which is all
 *	implicitly ones
 */
void smatrix_dsplud( smatrix_t *mat_in )
{
	long int i,j,k,m,idx,jdx,its,dx1,dx2;
	double xi,xk,tau,m1,m2;
	for(i=0;i<mat_in->ncols-1;i++)
	{
		smatrix_entry( mat_in, i, i, &xk );
		if( fabs(xk) < SMATRIX_ZERO_THRESHOLD )
		{
			printf( "smatrix_dsplud(): zero on the diagonal at (%d,%d) = %f\n", i, i, xk );
			exit(0);
		}
		for(j=i+1;j<mat_in->nrows;j++) // (j,j) is entry multiplied by unity
		{
			smatrix_entry( mat_in, j, i, &xi );
			tau = xi / xk;
			smatrix_insert( mat_in, j, i, tau );
			binary_closest( &i, mat_in->jvec[j], &mat_in->lvec[j], &jdx, &its );
			binary_closest( &i, mat_in->jvec[i], &mat_in->lvec[i], &idx, &its );
			for(k=jdx+1,m=idx+1;;)
			{
				m1 = 0.0; m2 = 0.0;
				dx1 = -1; dx2 = -1;
				if( k >= mat_in->lvec[j] && m >= mat_in->lvec[i] )
					break;
				if( k < mat_in->lvec[j] )
				{
					dx1 = mat_in->jvec[j][k];
					m1 = mat_in->mvec[j][k];
				}
				if( m < mat_in->lvec[i] )
				{
					dx2 = mat_in->jvec[i][m];
					m2 = mat_in->mvec[i][m];
				}
				if( dx1 != -1 && dx2 != -1 )
				{
					if( dx1 == dx2 )
					{
						mat_in->mvec[j][k] = m1 - tau * m2; // j,[k] entry guaranteed to exist
						++k,++m;
					}
					else
					{
						if( dx1 < dx2 )
							++k;
						else
						{
							smatrix_insert( mat_in, j, dx2, -1.0 * tau * m2 );
							++m;
						}
					}
				}
				else
				{
					if( dx1 > -1 ) // But dx2 is not
						++k;
					else
					{
						if( dx2 > -1 ) // But dx1 is not
						{
							smatrix_insert( mat_in, j, dx2, -1.0 * tau * m2 );
							++m,++k;
						}
						else // REACHING HERE IS A BUG
							printf( "smatrix_dsplud(): both indices i, j out of range; THIS IS BUG\n" );
					}
				}
			}
		}
	}
}

/**
 * This function takes as input the matrix output from smatrix_dsplud, which is the lower and upper
 * triangular matrices resulting from the decomposition.
 */
void smatrix_dspfbs_slow( smatrix_t *mat_in, smatrix_t *rhs_in )
{
	long int i,j,k,idx,jdx,kdx,its,cnt;
	int ret;
	double xji,ljk,bki,dnm;
	assert( mat_in->ncols == rhs_in->nrows );

	cnt = 0;
	for(i=0;i<rhs_in->ncols;i++)
	{
		for(j=0;j<rhs_in->nrows;j++)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=0;k<j;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, k, &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, k, i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_insert( rhs_in, j, i, xji );
		}
	}

	for(i=0;i<rhs_in->ncols;i++)
	{
		for(j=rhs_in->nrows-1;j>=0;j--)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=j+1;k<mat_in->ncols;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, k, &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, k, i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_entry( mat_in, j, j, &dnm );
			smatrix_insert( rhs_in, j, i, xji / dnm );
		}
	}
}

void smatrix_dspfbs( smatrix_t *mat_in, smatrix_t *rhs_in )
{
	long int i,j,k,idx,its;
	double xji,ljk,bki,dnm;
	assert( mat_in->ncols == rhs_in->nrows );

	for(i=0;i<rhs_in->ncols;i++)
	{
		for(j=0;j<rhs_in->nrows;j++)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=0;k<mat_in->lvec[j]&&mat_in->jvec[j][k]<j;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, mat_in->jvec[j][k], &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, mat_in->jvec[j][k], i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_insert( rhs_in, j, i, xji );
		}
	}

	for(i=0;i<rhs_in->ncols;i++)
	{
		for(j=rhs_in->nrows-1;j>=0;j--)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=j+1;k<mat_in->ncols;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, k, &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, k, i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_entry( mat_in, j, j, &dnm );
			smatrix_insert( rhs_in, j, i, xji / dnm );
		}
	}
}

void smatrix_dspfbs_size( smatrix_t *mat_in, smatrix_t *rhs_in, long int pos_in, long int size_in )
{
	long int i,j,k;
	double xji,ljk,bki,dnm;

	assert( mat_in->ncols == rhs_in->nrows );
	assert( pos_in < rhs_in->ncols && pos_in + size_in <= rhs_in->ncols );

	for(i=pos_in;i<pos_in+size_in;i++)
	{
		for(j=0;j<rhs_in->nrows;j++)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=0;k<j;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, k, &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, k, i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_insert( rhs_in, j, i, xji );
		}
	}

	for(i=pos_in;i<pos_in+size_in;i++)
	{
		for(j=rhs_in->nrows-1;j>=0;j--)
		{
			smatrix_entry( rhs_in, j, i, &xji );
			for(k=j+1;k<mat_in->ncols;k++)
			{
				/* l(j,k) */
				smatrix_entry( mat_in, j, k, &ljk );

				/* b(k,i) */
				smatrix_entry( rhs_in, k, i, &bki );

				/* Calculate x(j,i) */
				xji = xji - ljk * bki;
			}
			smatrix_entry( mat_in, j, j, &dnm );
			smatrix_insert( rhs_in, j, i, xji / dnm );
		}
	}
}

void smatrix_dspfbs_dvector( smatrix_t *mat_in, dvector_t *rhs_in )
{
	long int i,j,k,idx,jdx,kdx,its;
	int ret;
	double xji,ljk,bki,dnm;
	assert( mat_in->ncols == rhs_in->length );

	for(j=0;j<rhs_in->length;j++)
	{
// 		smatrix_entry( rhs_in, j, i, &xji );
		xji = rhs_in->start[j*rhs_in->rsize];
		for(k=0;k<j;k++)
		{
			/* l(j,k) */
			smatrix_entry( mat_in, j, k, &ljk );

			/* Calculate x(j,i) */
			xji = xji - ljk * rhs_in->start[k*rhs_in->rsize];
		}
		rhs_in->start[j*rhs_in->rsize] = xji;
	}

	for(j=rhs_in->length-1;j>=0;j--)
	{
		xji = rhs_in->start[j*rhs_in->rsize];
		for(k=j+1;k<mat_in->ncols;k++)
		{
			/* l(j,k) */
			smatrix_entry( mat_in, j, k, &ljk );

			/* Calculate x(j,i) */
			xji = xji - ljk * rhs_in->start[k*rhs_in->rsize];
		}
		smatrix_entry( mat_in, j, j, &dnm );
		rhs_in->start[j*rhs_in->rsize] = xji / dnm;
	}
}

void smatrix_darkrd_sparse_inv_prod( int n_in, smatrix_t *lud_in, smatrix_t *rhm_in, dvector_t *b_in, smatrix_t *Q_in, smatrix_t *H_in )
{
	int i,j,k;
	double sum,tmp;

	/* Copy b_in into first column of Q_in */
	sum = 0.0;
	for(i=0;i<b_in->length;i++)
		sum += b_in->start[i] * b_in->start[i];
	sum = sqrt( sum );
	for(i=0;i<b_in->length;i++)
		smatrix_insert( Q_in, i, 0, b_in->start[i] / sum );

	/* Subtract one from n_in to make basis of correct size */
	--n_in;

	/* Outer loop - Matrix-vector multiplication in Q */
	for(i=0;i<n_in;i++)
	{
		smatrix_dspmv_sparse( rhm_in, Q_in, i, Q_in, i + 1 );
		smatrix_dspfbs_size( lud_in, Q_in, i + 1, 1 );
		for(j=0;j<=i;j++)
		{
			tmp = smatrix_col_dot_prod( Q_in, j, Q_in, i + 1 );
			smatrix_lin_com( -1.0 * tmp, Q_in, j, 1.0, Q_in, i + 1, Q_in, i + 1 );
			smatrix_insert( H_in, j, i, tmp );
		}
		sum = smatrix_col_dot_prod( Q_in, i + 1, Q_in, i + 1 );
		sum = sqrt( sum );
		smatrix_insert( H_in, j, i, sum );
		for(k=0;k<Q_in->nrows;k++)
		{
			smatrix_entry( Q_in, k, i + 1, &tmp );
			smatrix_insert( Q_in, k, i + 1, tmp / sum );
		}
	}

	/* Create final column in H */
	smatrix_dspmv_sparse_dvector( rhm_in, Q_in, n_in, b_in );
	smatrix_dspfbs_dvector( lud_in, b_in );
	for(i=0;i<=n_in;i++)
	{
		tmp = smatrix_col_dot_prod_dvector( Q_in, i, b_in );
		smatrix_insert( H_in, i, n_in, tmp );
	}
}

