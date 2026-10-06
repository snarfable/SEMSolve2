#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <basis.h>
#include <mesh.h>
#include <ltgmap.h>
#include <semerror.h>
#include <dlinalg.h>
#include <slinalg.h>

#define RAND_RADIUS 2.0
#define INVERSE_MIN_RESIDUAL 1.0e-5
#define INVERSE_ITER_MAX 20000

void print_usage()
{
	printf( "Usage: semsolve [-n dim] [-d ideg] [-q qdeg] [-k krylov] [input mesh]\n" );
}

int main( int argc, char **argv )
{
	/* General garbage */
	int *ptr,yn,opt,krydim;
	long int i,j,k,m,d,q,ret,mdim,pdim,dim,row,col;
	double dbl;
	char tst[64];
	FILE *fp;

	if( argc < 2 )
	{
		print_usage();
		exit(0);
	}

	d = 1, q = 2, krydim = 100;

	while( ( opt = getopt( argc, argv, "n:d:q:k:" ) ) != -1 )
	{
		switch(opt)
		{
			case 'n':
				dim = atol( optarg );
				break;
			case 'd':
				d = atol( optarg );
				break;
			case 'q':
				q = atol( optarg );
				break;
			case 'k':
				krydim = atoi( optarg );
				break;
			default:
				print_usage();
				break;
		}
	}

	/* Solver and subspace projection variables */
	double evs[krydim];

	/* Mesh and basis variables */
	ltgmap_t map;
	basis_t basis;
	mesh_t mesh;
	potential_t pot;

	/* Shit for inverse iteration */
	dvector_t bv,tv;

	/* Linear algebra data structures */
	smatrix_t mat1,mat2,qm,hm;

	/* Initialize the basis */
	basis_init_std( &basis, 3, d, q );

	/* Load the mesh */
	mesh_init( &mesh, 3 );
	mesh_load( &mesh, argv[optind], 1 );
	mesh_gen_node_simp_map( &mesh );

	/* Load the ions */
	int chgs[1] = { 1 };
	double ions[3] = { 5.0, 5.0, 5.0 };

	/* Put the ion data into a potential_t structure */
	potential_init( &pot, 3, 1 );
	bcopy( chgs, pot.charges, pot.nions * sizeof(int) );
	bcopy( ions, pot.coords, pot.nions * pot.dim * sizeof(double) );

	/* Initialize the local-to-global map */
	ltgmap_init( &map, &basis, &mesh, &pot );

	/* Output to make sure everything is pretty */
	printf( "nglobal = %d\n", map.nglobal );
	printf( "bglobal = %d\n", map.bglobal );

	/* Set the system dimension */
	mdim = map.bglobal;
	pdim = map.basis->pdim;
	dim = map.basis->dim;

	bv.start = (double*) malloc( mdim * sizeof(double) );
	bv.rsize = 1;
	bv.length = mdim;

	srand((unsigned)time(0));
	for(i=0;i<mdim;i++)
		bv.start[i*bv.rsize] = RAND_RADIUS * (rand() > RAND_MAX / 2 ? -1.0 : 1.0) * (double) rand() / (double) RAND_MAX;

	/* Now make the matrices */
	smatrix_init( &mat1, mdim, mdim );
	smatrix_init( &mat2, mdim, mdim );
	smatrix_init( &qm, mdim, krydim );
	smatrix_init( &hm, krydim, krydim );

	for(i=0;i<map.mesh->nsimp;i++)
	{
		for(j=0;j<map.basis->pdim;j++)
			for(k=j;k<map.basis->pdim;k++)
		{
			row = map.map[i*pdim+j];
			col = map.map[i*pdim+k];
			if( map.bexcld[row] == -1 || map.bexcld[col] == -1 )
				continue;
			else
			{
				row = map.bexcld[row];
				col = map.bexcld[col];
			}
			dbl = ltgmap_inner_prod( &map, i, j, k );
			smatrix_entry_add( &mat1, row, col, dbl );
			if(row!=col)
				smatrix_entry_add( &mat1, col, row, dbl );
			dbl = 0.5 * ltgmap_kinetic_inner_prod( &map, i, j, k ) + ltgmap_potential_inner_prod( &map, i, j, k );
			smatrix_entry_add( &mat2, row, col, dbl );
			if(row!=col)
				smatrix_entry_add( &mat2, col, row, dbl );
		}
	}

	/* Outputting */
	smatrix_print_debug( &mat2 );

	/* Solve for relevant eigenvalues by subspace projection */
	smatrix_dsplud( &mat1 );
	smatrix_darkrd_sparse_inv_prod( krydim, &mat1, &mat2, &bv, &qm, &hm );

	/* Then solve */
	smatrix_eigen_solve( &hm );

	/* Output eigenvalues */
	for(i=0;i<krydim;i++)
		smatrix_entry( &hm, i, i, evs + i );

	/* Sort eigenvalues */
	bubble_sort_double_long( krydim, evs );

	/* Save the shit */
	fp = fopen( "eigs.out", "w" );
	for(i=0;i<krydim;i++)
		fprintf( fp, "%f ", evs[i] );
	fclose( fp );

	/* Now ask whether to continue... */
	for(;;)
	{
		printf( "Continue to inverse iteration? (y/n): " );
		scanf( "%s", tst );
		if( tst[0] == 'n' || tst[0] == 'N' )
			exit(0);
		if( tst[0] == 'y' || tst[0] == 'Y' )
			break;
	}

	/* Now start inverse iteration */
	printf( "Inverting the mass matrix: " );
	fflush( stdout );
	smatrix_dspfbs( &mat1, &mat2 ); /* Now mat2 contains the matrix */
	printf( "Done.\n" );
	fflush( stdout );
	printf( "Shifting: " );
	fflush( stdout );
	smatrix_shift( &mat2, evs[0] ); /* Shift by the bottom eigenvalue */
	printf( "Done.\n" );
	fflush( stdout );

	printf( "LU Decomposing total system: " );
	fflush( stdout );
	smatrix_dsplud( &mat2 );
	printf( "Done.\n" );
	fflush( stdout );

	printf( "Inverse iteration: " );
	fflush( stdout );
	bv.start = (double*) malloc( mdim * sizeof(double) );
	bv.rsize = 1;
	bv.length = mdim;
	tv.start = (double*) malloc( mdim * sizeof(double) );
	tv.rsize = 1;
	tv.length = mdim;
	dbl = 0.0;
	srand((unsigned)time(0));
	for(i=0;i<mdim;i++)
	{
		bv.start[i*bv.rsize] = (rand() > RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
		dbl += bv.start[i*bv.rsize] * bv.start[i*bv.rsize];
	}
	dbl = sqrt( dbl );
	for(i=0;i<mdim;i++)
		bv.start[i*bv.rsize] /= dbl;
	for(i=0;i<INVERSE_ITER_MAX;i++)
	{
		for(j=0;j<bv.length;j++) // Copy into tv first
			tv.start[j*tv.rsize] = bv.start[j*bv.rsize];
		smatrix_dspfbs_dvector( &mat2, &bv );
		dbl = 0.0;
		for(j=0;j<mdim;j++)
			dbl += bv.start[j*bv.rsize] * bv.start[j*bv.rsize];
		dbl = sqrt( dbl );
		for(j=0;j<mdim;j++)
			bv.start[j*bv.rsize] /= dbl; // Normalize to unit vector
		dbl = 0.0;
		for(j=0;j<mdim;j++)
			dbl += powf( tv.start[j*tv.rsize] - bv.start[j*bv.rsize], 2.0 );
		dbl = sqrt( dbl );
		printf( "r = %f\n", dbl ); // Print out the residual if verbose
		fflush( stdout );
		if( dbl < INVERSE_MIN_RESIDUAL )
			break;
	}
	printf( "Done.\n" );
	fflush( stdout );

	printf( "Final vector: " );
	for(i=0;i<mdim;i++)
		printf( "%f ", bv.start[i*bv.rsize] );
	printf( "\n" );
	fflush( stdout );

	/* Start the plotting routine */
	printf( "Plotting eigenvector: " );
	fflush( stdout );

	mesh_write( &mesh, "plot.msh", 0 );

	fp = fopen( "plot.msh", "a" ); /* Append the node data on the end */
	if( fp == NULL )
	{
		printf( "Could not open output for writing\n" );
		exit(0);
	}
	fprintf( fp, "$NodeData\n" );
	fprintf( fp, "1\n" );		// 1 string tag
	fprintf( fp, "\"PLOT1\"\n" );	// View name
	fprintf( fp, "1\n" );		// 1 scalar tag
	fprintf( fp, "0.0\n" );		// The time
	fprintf( fp, "3\n" );		// 3 integer tags
	fprintf( fp, "0\n" );		// The time step
	fprintf( fp, "1\n" );		// 1-component vector field
	fprintf( fp, "%d\n", mesh.nvert );
	fflush( fp );

	for(i=0;i<mesh.nvert;i++)
	{
		k = mesh.nsmap[i][0]; /* Index of the first simplex containing vertex i; there will always be at least one */
		m = simple_find_long( i, dim + 1, mesh.simpx + k * ( dim + 1 ) );
		if( m == -1 )
		{
			printf( "The node-simplex map is inconsistent. Exiting.\n" );
			exit(0);
		}
		k = map.map[k*basis.pdim+m];
		m = map.bexcld[k];
		if( m == -1 )
			fprintf( fp, "%d 0.0\n", i + 1 );
		else
			fprintf( fp, "%d %f\n", i + 1, bv.start[m*bv.rsize] );
	}
	fprintf( fp, "$EndNodeData\n" );
	fclose( fp );

	printf( "Done.\n" );
	fflush( stdout );

	smatrix_free( &mat1 );
	smatrix_free( &mat2 );

	basis_free( &basis );
	mesh_free( &mesh );
	ltgmap_free( &map );

	free( bv.start );

	return 0;
}
