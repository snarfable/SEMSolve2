#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <unistd.h>

int main( int argc, char **argv )
{
	if( argc < 5 )
	{
		printf( "%s: Too few arguments\n", argv[0] );
		exit(0);
	}

	int i,j;
	int vdim = atoi( argv[1] );
	int nvec = atoi( argv[2] );
	double rrad = atof( argv[3] );
	double tmp[nvec*vdim];
	double sum;
	long int n;
	FILE *fp;

	fp = fopen( argv[4], "w" );
	if( fp == NULL )
	{
		printf( "%s: Could not open file `%s' for writing\n", argv[0], argv[4] );
		exit(0);
	}
	srand((unsigned)time(0));
	n = 0;
	while( n < nvec )
	{
		sum = 0.0;
		for(j=0;j<vdim;j++)
		{
			tmp[n*vdim+j] = (rand()<RAND_MAX/2?1.0:-1.0) * rrad * (double) rand() / (double) RAND_MAX;
			sum += tmp[n*vdim+j] * tmp[n*vdim+j];
		}
		sum = sqrt( sum );
		if( sum <  rrad )
			++n;
	}

	fprintf( fp, "%d\n%d\n", vdim, nvec );
	for(i=0;i<nvec;i++)
	{
		for(j=0;j<vdim;j++)
			fprintf( fp, "%f ", tmp[i*vdim+j] );
		fprintf( fp, "\n" );
	}

	fclose( fp );
	return 0;
}
