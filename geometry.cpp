#include <geometry.h>

double length( int dim_in, double *point_in )
{
	int i;
	double tot = 0.0;
	for(i=0;i<dim_in;i++)
		tot += point_in[i] * point_in[i];
	return sqrt( tot );
}

double metric( int dim_in, double *vtx1_in, double *vtx2_in )
{
	int i;
	double tot = 0.0;
	for(i=0;i<dim_in;i++)
		tot += powf( vtx1_in[i] - vtx2_in[i], 2.0 );
	return sqrt( tot );
}

/**
 * WARNING: This function changes the inputs!!!
 */
void normalize( int dim_in, double *vtx_in )
{
	int i;
	double len;
	len = 0.0;
	for(i=0;i<dim_in;i++)
		len += vtx_in[i] * vtx_in[i];
	len = sqrt( len );
	for(i=0;i<dim_in;i++)
		vtx_in[i] /= len;
}

void vec_sub( int dim_in, double *vtx1_in, double *vtx2_in )
{
	int i;
	for(i=0;i<dim_in;i++)
		vtx2_in[i] -= vtx1_in[i];
}

/**
 * WARNING: This function changes the inputs!!!
 */
double angle( int dim_in, double *u1_in, double *u2_in, double *v1_in, double *v2_in )
{
	int i;
	double dot;
	for(i=0;i<dim_in;i++)
		u2_in[i] -= u1_in[i], v2_in[i] -= v1_in[i];
	normalize( dim_in, u2_in );
	normalize( dim_in, v2_in );
	dot = dotp( dim_in, u2_in, v2_in );
	return acos( dot );
}

/**
 * This function takes vector @e vec_in @e - @e orig_in and projects
 * it onto @e norm_in @e - @e orig_in , and stores the output in
 * @e norm_in .
 */
void project( int dim_in, double *orig_in, double *norm_in, double *vec_in )
{
	int i;
	double dp;
	vec_sub( dim_in, orig_in, norm_in );
	vec_sub( dim_in, orig_in, vec_in );
	normalize( dim_in, norm_in );
	dp = dotp( dim_in, norm_in, vec_in );
	for(i=0;i<dim_in;i++)
		norm_in[i] *= dp;
}

/**
 * Notice that the origin is assumed to be the first vector of the manifold. It must be used
 * in order for everything to be consistent.
 */
void subspace_project( int dim_in, int mdim_in, double *mnfd_in, double *point_in, double *proj_out )
{
	int i,j;
	double dot,len;

	/* Keep the running total */
	for(i=0;i<dim_in;i++)
		proj_out[i] = 0.0;

	/* Want point_in relative to mnfd_in[0] */
	vec_sub( dim_in, mnfd_in, point_in );
	for(i=1;i<mdim_in;i++)
		vec_sub( dim_in, mnfd_in, mnfd_in + i * dim_in );

	/* Now project onto each generalized surface vertex */
	for(i=1;i<mdim_in;i++)
	{
		len = length( dim_in, mnfd_in + i * dim_in );
		dot = dotp( dim_in, point_in, mnfd_in + i * dim_in ) / len;
		for(j=0;j<dim_in;j++)
			proj_out[j] += dot * mnfd_in[i*dim_in+j] / len;
	}
	for(i=0;i<dim_in;i++)
		proj_out[i] += mnfd_in[i];
}

/**
 * This function is dimension-generalized. @e angles_in must contain @e dim_in @e - @e 2
 * rotation parameters for arbitrary dimension. Of course, with @e dim_in @e = @e 3 , only
 * a single angle is needed.
 */
void rotate_axis( int dim_in, int mdim_in, double *mnfd_in, double *point_in, double angle_in, double *point_out )
{
	int i,j;
	double n[dim_in],p[dim_in],q[dim_in],dot,len,x,y;

	/* Preserve original vectors */
	bcopy( point_in, p, dim_in * sizeof(double) );

	/* Calculate p = r - r proj n */
	subspace_project( dim_in, mdim_in, mnfd_in, point_in, n );
	vec_sub( dim_in, n, p ); /* p is now a vector pointing away from the manifold toward point_in */
	len = length( dim_in, p ); /* len is the radius of point_in from the manifold */
	normalize( dim_in, p ); /* Now point_in = n + len * p */

	/* Now do Gramm-Schmidt to orthogonalize relative to mnfd_in and p */
	for(i=0;i<dim_in;i++)
		q[i] = 1.0; /* Does not matter what vector is used */
	for(i=1;i<mdim_in;i++)
	{
		normalize( dim_in, mnfd_in + i * dim_in );
		dot = dotp( dim_in, mnfd_in + i * dim_in, q );
		for(j=0;j<dim_in;j++)
			q[j] -= dot * mnfd_in[i*dim_in+j];
	}

	/* FIXME: Need to standardize the DIRECTION, i.e. + or -, of q to make input angles correct	*/
	/* FIXME: Want vectors listed in right-handed order relative to the manifold			*/
	dot = dotp( dim_in, p, q );
	for(i=0;i<dim_in;i++)
		q[i] -= dot * p[i];
	assert( length( dim_in, q ) > 0.0 );
	normalize( dim_in, q );

	/* Now p (adjacent) and q (opposite) are the vectors which determine the rotation */
	for(i=0;i<dim_in;i++)
		point_out[i] = 0.0;
	x = cos( angle_in );
	y = sin( angle_in );
	for(i=0;i<dim_in;i++)
		point_out[i] += len * ( x * p[i] - y * q[i] ) + n[i];
}
