#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mesh.h>

void print_usage()
{
	printf( "Usage: convert_mesh [mesh input] [in format #] [mesh output] [out format]\n" );
}

int main( int argc, char **argv )
{
	int i;
	mesh_t mesh;

	if( argc < 5 )
	{
		print_usage();
		exit(0);
	}

// 	mesh_load( &mesh, argv[1],  );

	return 0;
}
