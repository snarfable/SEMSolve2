#include <ltgmap.h>
#include <stdio.h>

void print_usage()
{
	printf( "Usage: mk_mesh [dimension] [nodes file] [output file] [format]\n" );
}

int main( int argc, char **argv )
{
	mesh_t mesh;

	if( argc < 5 )
	{
		printf( "%s: Too few arguments\n", argv[0] );
		print_usage();
		return 0;
	}

	/* Initialize shit */
	mesh_init( &mesh, atoi( argv[1] ) );
	mesh_load_nodes( &mesh, argv[2], 0 );
	mesh_gen_convex_hull( &mesh );
// 	mesh_gen_del_wrap( &mesh );

	/* Write the mesh to file */
// 	mesh_write( &mesh, argv[3], atoi( argv[4] ) );

	/* Free the shit */
	mesh_free( &mesh );

	return 0;
}
