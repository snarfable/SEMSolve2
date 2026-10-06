/**
 * This file contains the implementations of the utility functions for
 * the dense linear algebra header dlinalg.h
 * \file slinalg_util.cpp
 * \author Ryan Daly
 * \em rpdaly@ucsd.edu
 * \date 03/09/2009
 */

#include <dlinalg.h>

void dmatrix_print( dmatrix_t *mat_in )
{
	long int i,j;
	for(i=0;i<mat_in->nrows;i++)
	{
		for(j=0;j<mat_in->ncols;j++)
			printf( "%f ", mat_in->start[i*mat_in->rsize+j] );
		printf( "\n" );
	}
}
