#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <err.h>
#include <combinadic.h>
#include <ltgmap.h>
#include <mesh.h>
#include <geometry.h>

#define SMOOTH_MESH_MAX_ITS 3

void print_usage()
{
	printf( "Usage: smooth_mesh [dim] [input mesh] [output mesh] [fmt] [num its]\n" );
}

void error_wrapper( int ret_in )
{
	if( ret_in == -1 )
	{
		fprintf( stderr, "Error, exiting\n" );
		exit(0);
	}
}

int main( int argc, char **argv )
{
	int i,j,k,m,n,p,c,dim,bad,*bnd;
	long int f;
	double tmp,*ptr,*loc,*avg,*ivols = NULL,step = 0.001;

	mesh_t mesh;
	ltgmap_t map;
	list_t list;

	if( argc < 6 )
	{
		printf( "%s: Too few arguments\n", argv[0] );
		print_usage();
		exit(0);
	}

	/* Number of iterations to take */
	int nits = atoi( argv[5] );

	/* Read shit from command line */
	dim = atoi( argv[1] );
	avg = (double*) malloc( dim * sizeof(double) );
	loc = (double*) malloc( dim * sizeof(double) );

	/* Load the mesh */
	mesh_init( &mesh, dim );
	mesh_load( &mesh, argv[2], 1 );

	/* Initialize temp storage */
	double nodes[dim*mesh.nvert];

	/* Need to load ltgmap_t::nsmap since it is needed by ltgmap_gen_planes_set() */
	mesh_gen_node_simp_map( &mesh );

	/* Only generating planes set, i.e. plane, member, psize, and mbsize */
	ltgmap_init_mesh_only( &map, &mesh );

	/* Detect boundary nodes */
	bnd = (int*) malloc( mesh.nvert * sizeof(int) );
	if( bnd == NULL )
	{
		fprintf( stderr, "Could not allocate memory\n" );
		exit(0);
	}

	/* Initialize all entries to 0 (NOT on boundary) */
	for(i=0;i<mesh.nvert;i++)
		bnd[i] = 0;

	/* Iterate all planes in the mesh */
	for(i=0;i<map.nplane;i++)
	{
		/* Look for simplices of order dim */
		if( map.psize[i] != dim )
			continue;

		/* Make sure only one simplex contains this face */
		if( map.mbsize[i] > 1 )
			continue;

		/* Otherwise set nodes contained therein to boundary */
		for(j=0;j<dim;j++)
			if( bnd[map.plane[i][j]] == 0 )
				bnd[map.plane[i][j]] = 1;
	}

	for(c=0;c<nits;c++)
	{
		/* For each vertex */
		for(i=0;i<mesh.nvert;i++)
		{
			/* Do NOT update boundary vertices */
			if( bnd[i] == 1 )
				bcopy( mesh.vertx + i * dim, nodes + i * dim, dim * sizeof(double) ); /* Copy the current location into nodes */

			/* Zero the centroid */
			for(j=0;j<dim;j++)
				avg[j] = 0.0;

			/* Now measure initial volumes */
			ptr = (double*) realloc( ivols, mesh.nssize[i] * sizeof(double) );
			if( ptr == NULL )
			{
				fprintf( stderr, "%d - %d - ", mesh.nssize[i], i );
				fprintf( stderr, "Memory error\n" );
// 				exit(0);
			}
			ivols = ptr;
			for(j=0;j<mesh.nssize[i];j++)
				ivols[j] = ltgmap_simplex_volume_signed( &map, mesh.nsmap[i][j] );

			/* For each simplex containing i */
			error_wrapper( list_init( &list, sizeof(long int) ) );
			for(j=0;j<mesh.nssize[i];j++)
			{
				/* Select a simplex */
				m = mesh.nsmap[i][j];

				/* Project vertex i onto the face and measure its distance */
				for(k=0;k<dim+1;k++)
				{
					n = mesh.simpx[m*(dim+1)+k];
					f = simple_find_long( n, list.size, (long int*) list.start );
					if( f == -1 ) /* If not found in list_in */
						error_wrapper( list_append( &list, (void*) &n ) );
					else
						continue;
					for(p=0;p<dim;p++)
						avg[p] += mesh.vertx[n*dim+p];
				}
			}

			/* Divide by number of unique vertices */
			for(j=0;j<dim;j++)
				avg[j] /= (double) list.size; /* Now avg is the centroid */

			/* Free the list */
			error_wrapper( list_free( &list ) );

			/* Calculate the gradient; put it in avg */
			for(j=0;j<dim;j++)
				avg[j] -= mesh.vertx[i*dim+j];

			/* Copy nodes(i) into loc before updating in case this is a bad move */
			bcopy( mesh.vertx + i * dim, loc, dim * sizeof(double) );

			/* Offset the node in direction of avg using step */
			for(j=0;j<dim;j++)
				nodes[i*dim+j] = mesh.vertx[i*dim+j] + step * avg[j];

			/* Swap ltgmap mesh.vertx for nodes temporarily for volume calculation */
			ptr = mesh.vertx;
			mesh.vertx = nodes;

			/* Now measure final volumes and compare to initials */
			bad = 0;
			for(j=0;j<mesh.nssize[i];j++)
			{
				tmp = ltgmap_simplex_volume_signed( &map, mesh.nsmap[i][j] );
				if( ( tmp < 0.0 && ivols[j] > 0.0 ) || ( tmp > 0.0 && ivols[j] < 0.0 ) )
				{
					bad = 1;
					break;
				}
			}

			if( bad == 1 )
			{
// 				fprintf( stdout, "Bad move; volume sign change in element %d\n", mesh.nsmap[i][j] );
				bcopy( loc, nodes + i * dim, dim * sizeof(double) );
			}

			/* Revert to original vertx */
			mesh.vertx = ptr;
		}
		bcopy( nodes, mesh.vertx, dim * mesh.nvert * sizeof(double) );
	}
	free( avg );
	free( loc );

	/* Write the mesh */
	mesh_write( &mesh, argv[3], atoi( argv[4] ) );

	/* Bye bye */
	return 0;
}
