#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlinalg.h>
#include <basis.h>
#include <mesh.h>
#include <ltgmap.h>
#include <geometry.h>

#define HARD_BOUNDARY

void print_usage()
{
	printf( "Usage: smooth_dist [dim] [node file in] [node file out] [num iterations]\n" );
}

int main( int argc, char **argv )
{
	if( argc < 5 )
	{
		print_usage();
		exit(0);
	}

	int i,j,k,m,dim = atoi( argv[1] );
	double frc[dim],tfc[dim],ctr[] = { 5.0, 5.0, 5.0 };
	mesh_t mesh;
	FILE *fp;

	mesh_init( &mesh, (long) dim );
	if( mesh_load_nodes( &mesh, argv[2], 0 ) != 0 )
	{
		fprintf( stderr, "%s: Cannot load nodes from file %s\n", argv[0], argv[2] );
		return 0;
	}

	for(i=0;i<atoi( argv[4] );i++)
	{
		
	}

	fp = fopen( argv[3], "w" );
	if( fp == NULL )
	{
		fprintf( stderr, "%s: Cannot save nodes to file %s\n", argv[0], argv[3] );
		return 0;
	}

	/* Output some garbage */
	fprintf( fp, "%d\n%d\n", mesh.dim, mesh.nvert );
	for(i=0;i<mesh.nvert;i++)
	{
		for(j=0;j<mesh.dim;j++)
			fprintf( fp, "%f ", mesh.vertx[i*mesh.dim+j] );
		fprintf( fp, "\n" );
	}
	fclose( fp );

	return 0;
}
