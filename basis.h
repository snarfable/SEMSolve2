#ifndef SEMSOLVE2_BASIS_H
#define SEMSOLVE2_BASIS_H

#define SEMSOLVE2_BASIS_USE_LAPACK

/**
 * This is the header file which defines the basis_t type/class. Version 2 of this
 * is more simplistic than the class-based older version.
 * \file basis.h
 * \author Ryan Daly
 * \em rpdaly@ucsd.edu
 * \date 03/04/2009
 */

#include <combinadic.h>
#include <math.h>

#include <lapack.h>
#include <dlinalg.h>

typedef struct
{
	long int dim;
	long int deg;
	long int pdim;
	double	*coeffs;
	double	*nodes;
	long int qdeg;
	long int qpdim;
	double	*qcoeffs;
	double	*qnodes;
} basis_t;

/*
 * Creation/destruction
 */
void	 basis_init( basis_t *, long int, long int, long int );
void	 basis_free( basis_t * );
void	 basis_init_std( basis_t *, long int, long int, long int );
void	 basis_load_nodes( basis_t *, char * );
void	 basis_load_std_nodes( basis_t * );						///< Must have called basis_init first
void	 basis_load_quad_nodes( basis_t *, char * );					///< Must have initialized so dim, deg, and qdeg have values
void	 basis_load_std_quad_nodes( basis_t * );

/*
 * Output functions for additional irritation
 */
void basis_print_interp_summary( basis_t * );
void basis_print_quad_summary( basis_t * );

/*
 * Functions for general annoyance and bookkeeping
 */
long int basis_term_index( long int, long int, long int * );
int	 basis_term_exponent( long int, long int, long int, long int * );
long int basis_poly_dim( long int, long int );

/*
 * Mathematical functions
 */

/* General */
int	 basis_gen_coeffs( basis_t * );							///< This function is called when the basis is constructed
int	 basis_gen_quad_coeffs( basis_t * );						///< This function is called when the basis is constructed
double	 basis_term_integral_simplex( long int, long int * );				///< A very very nice function for integrating single polynomial terms over simplices
double	 basis_plane_test( long int, long int, double *, long int *, double * );	///< Used to detect whether an input node lies in a given plane in the simplex
int	 basis_func_gradient( basis_t *, long int, double *, double * );		///< basis, func index, input, output
void	 basis_barycentric( long int, double *, double *, double * );			///< Used to calculate the barycentric coordinates with a d-simplex as polygon

/* Evaluation of basis functions at a point */
double	 basis_eval_interp_func( basis_t *, long int, double * );			///< Evaluates an interpolation function at a given input point

#endif
