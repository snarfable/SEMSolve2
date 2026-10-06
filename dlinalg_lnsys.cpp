#include <dlinalg.h>

double dvector_dotp( dvector_t *v1_in, dvector_t *v2_in )
{
	assert( v2_in->length == v1_in->length );

	long int i;
	double ret = 0.0;

	for(i=0;i<v1_in->length;i++)
		ret += v1_in->start[i*v1_in->rsize] * v2_in->start[i*v2_in->rsize];

	return ret;
}

double dotp( int dim_in, double *v1_in, double *v2_in )
{
	int i;
	double tot = 0.0;

	for(i=0;i<dim_in;i++)
		tot += v1_in[i] * v2_in[i];

	return tot;
}

/**
 * This function outputs a pivot vector which identifies the rows
 * such that there are no zeroes on the diagonal of the matrix.
 */
void dmatrix_pivot( dmatrix_t *mat_in, int *piv_out )
{
	
}

void dmatrix_lu_dcmp( dmatrix_t *mat_in )
{
	long int i,j,k;
	long int size = mat_in->nrows < mat_in->ncols ? mat_in->nrows : mat_in->ncols;
	double tau;
	for(i=0;i<size;i++) // mat_in(i,i), size => nrows
	{
		for(j=i+1;j<mat_in->nrows;j++) // mat_in(i+1,i), mat_in(i+2,i), ...
		{
			tau = mat_in->start[j*mat_in->rsize+i] / mat_in->start[i*mat_in->rsize+i];
			for(k=i+1;k<mat_in->ncols;k++)
				mat_in->start[j*mat_in->rsize+k] = mat_in->start[j*mat_in->rsize+k]
									- tau * mat_in->start[i*mat_in->rsize+k];
			mat_in->start[j*mat_in->rsize+i] = tau;
		}
	}
}

void dmatrix_lu_dcmp_bl( dmatrix_t *mat_in, long int bls_in )
{
	
}

void dmatrix_fb_sub_vrhs( dmatrix_t *mat_in, dvector_t *rhs_in )
{
	long int i,j;
	register double sum;
	assert( rhs_in->length == mat_in->ncols && mat_in->nrows == mat_in->ncols );
	int size = mat_in->nrows;
	for(i=1;i<size;i++)
	{
		sum = rhs_in->start[i*rhs_in->rsize];
		for(j=0;j<i;j++)
			sum -= mat_in->start[i*mat_in->rsize+j] * rhs_in->start[j*rhs_in->rsize];
		rhs_in->start[i*rhs_in->rsize] = sum;
	}
	for(i=size-1;i>=0;i--)
	{
		sum = rhs_in->start[i*rhs_in->rsize];
		for(j=size-1;j>i;j--)
			sum -= mat_in->start[i*mat_in->rsize+j] * rhs_in->start[j*rhs_in->rsize];
		sum /= mat_in->start[i*mat_in->rsize+i];
		rhs_in->start[i*rhs_in->rsize] = sum;
	}
}

void dmatrix_fb_sub_mrhs( dmatrix_t *mat_in, dmatrix_t *rhs_in )
{
	
}
