#include <potential.h>

/* XML Dependencies */
#include <libxml2/libxml/parser.h>

int potential_init( potential_t *pot_in, int dim_in, int nion_in )
{
	int *tmp;
	double *ump;

	pot_in->dim = dim_in;
	pot_in->nions = nion_in;

	tmp = (int*) malloc( nion_in * sizeof(int) );
	if( tmp == NULL )
		return -1;
	pot_in->charges = tmp;

	ump = (double*) malloc( dim_in * nion_in * sizeof(double) );
	if( ump == NULL )
		return -2;
	pot_in->coords = ump;

	return 0;
}

int potential_free( potential_t *pot_in )
{
	free( pot_in->charges );
	free( pot_in->coords );
}

int potential_load( potential_t *pot_in, char *file_in )
{
	
}

double potential_energy( potential_t *pot_in, double *pos_in )
{
	long int i,j;
	double tmp,tot = 0.0;

	/* For each nucleus */
	for(i=0;i<pot_in->nions;i++)
	{
		tmp = 0.0;
		for(j=0;j<pot_in->dim;j++)
			tmp += powf( pos_in[j] - pot_in->coords[i*pot_in->dim+j], 2.0 );
		tot -= 1.0 / sqrt( tmp );
	}
// 	printf( "tot = %f\n", tot );
	return tot;
}
