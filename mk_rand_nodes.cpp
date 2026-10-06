#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#include <mesh.h>
#include <basis.h>
#include <ltgmap.h>
#include <slinalg.h>
#include <dlinalg.h>
#include <semerror.h>

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
	FILE *fp;

	fp = fopen( argv[4], "w" );
	if( fp == NULL )
	{
		printf( "%s: Could not open file `%s' for writing\n", argv[0], argv[4] );
		exit(0);
	}
	fprintf( fp, "%d\n%d\n", vdim, nvec );
	srand((unsigned)time(0));
	for(i=0;i<nvec;i++)
	{
		for(j=0;j<vdim;j++)
			fprintf( fp, "%f ", (rand()<RAND_MAX/2?1.0:-1.0) * rrad * (double) rand() / (double) RAND_MAX );
		fprintf( fp, "\n" );
	}

	fclose( fp );
	return 0;
}
