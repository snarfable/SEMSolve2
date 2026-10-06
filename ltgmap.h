#ifndef SEMSOLVE2_LTGMAP_H
#define SEMSOLVE2_LTGMAP_H

#include <assert.h>
#include <math.h>
#include <basis.h>
#include <mesh.h>
#include <integrals.h>
#include <potential.h>

#define LTGMAP_ZERO_THRESHOLD 1.0e-7

typedef struct
{
	basis_t     *basis;			///< The local element data
	mesh_t      *mesh;			///< Yep, the mesh associated with this map
	potential_t *potential;			///< Potential energy function for nuclei
	long int    nglobal;			///< The total number of unique global nodes
	long int    *map;			///< The very reason we are here
	long int    nplane;			///< The number of planes stored in psize and plane
	long int    *mbsize;			///< Number of simplices containing each plane
	long int    *psize;			///< List of plane sizes, which can be from 1 up to dim, for use in plane
	long int    **member;			///< List of simplex indices which share the plane with the same index
	long int    **plane;			///< List of all unique planes
	long int    **npmap;			///< Maps local interpolation nodes to local plane combinadic indices
	long int    *npsize;			///< The number of interpolation nodes in each plane, referring to npmap
	long int    **npeqvs;			///< Interfacially equivalent nodes
	long int    *eqsize;			///< The number of nodes in each
	long int    nnpeqv;			///< The number of plane pairs, totalled
	long int    nnpmap;			///< The number of entries in npmap and npsize
	long int    *bexcld;			///< Final global map excluding the boundary nodes
	long int    bglobal;			///< The total number of nodes on which calculation should be done, actually
	char        *nbndry;			///< Is 0 if node i is not on boundary and 1 if is
} ltgmap_t;

/*
 * Creating and destroying
 */
void ltgmap_init( ltgmap_t *, basis_t *, mesh_t *, potential_t * );
void ltgmap_init_mesh_only( ltgmap_t *, mesh_t * );
void ltgmap_free( ltgmap_t * );

/*
 * For debugging
 */
void ltgmap_print_mesh_planes( ltgmap_t * );
void ltgmap_print_local_node_map( ltgmap_t * );
void ltgmap_print_face_map( ltgmap_t * );
void ltgmap_print_map( ltgmap_t * );
void ltgmap_print_boundary( ltgmap_t * );

/*
 * Math stuff
 */
void   ltgmap_affine( ltgmap_t *, long int, double *, double * );
double ltgmap_simplex_volume( ltgmap_t *, long int );
double ltgmap_simplex_volume_signed( ltgmap_t *, long int );
double ltgmap_inner_prod( ltgmap_t *, long int, long int, long int );
double ltgmap_potential_inner_prod( ltgmap_t *, long int, long int, long int );
double ltgmap_kinetic_inner_prod( ltgmap_t *, long int, long int, long int );

/*
 * Important functions
 */
int      ltgmap_gen_planes_set( ltgmap_t * );
int      ltgmap_gen_local_node_plane_map( ltgmap_t * );	///< Fills in the details for npmap
int      ltgmap_gen_local_plane_equivs( ltgmap_t * );	///< For gen'ing equivalencies b/w local face nodes
long int ltgmap_face_pair_index( long int, long int, long int, long int, long int, long int );
int      ltgmap_gen_ltgmap( ltgmap_t * );		///< Fills in the details for ltgmap_t::map
int	 ltgmap_boundary_check( ltgmap_t *, double * );	///< Takes a node anywhere in a mesh and detects whether or not it is on the boundary; this is EXPENSIVE so don't use it unless necessary
int	 ltgmap_detect_boundary( ltgmap_t * );		///< Fills in the nbndry array
int	 ltgmap_exclude_boundary( ltgmap_t * );		///< Fills in the bexcld array
int      ltgmap_simplex_local_indices( ltgmap_t *, long int, long int, long int *, long int * );

#endif
