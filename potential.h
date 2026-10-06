#ifndef SEMSOLVE2_POTENTIAL_H
#define SEMSOLVE2_POTENTIAL_H

#include <stdlib.h>
#include <math.h>

typedef struct
{
	int	dim;
	int	nions;
	int	*charges;
	double	*coords;
} potential_t;

/* Construction/destruction/loading */
int    potential_init( potential_t *, int, int );
int    potential_free( potential_t * );
int    potential_load( potential_t *, char * );

/* Specific energy functions */
double potential_energy( potential_t *, double * );

#endif
