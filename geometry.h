#ifndef SEMSOLVE2_GEOMETRY_H
#define SEMSOLVE2_GEOMETRY_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <dlinalg.h>

double length( int, double * );							///< Calculate the magnitude of a vector
double metric( int, double *, double * );					///< Calculate the Euclidean distance b/w two points
void   normalize( int, double * );						///< Normalizes a vector
void   vec_sub( int, double *, double * );					///< Subtract vector two from vector one
double angle( int, double *, double *, double *, double * );			///< Calculate the angle between two vectors
void   project( int, double *, double *, double * );				///< Project a vector onto another vector
void   subspace_project( int, int, double *, double *, double * );		///< Project a point onto a manifold of arbitrary dim
void   rotate_axis( int, int, double *, double *, double, double * );		///< Rotate a point about an arbitrary axis

#endif
