#ifndef SEMSOLVE2_MESH_H
#define SEMSOLVE2_MESH_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <combinadic.h>
#include <mesh.h>
#include <dlinalg.h>

#define SEMSOLVE2_USE_LAPACK 0
#define SEMSOLVE2_ALLOC_DATA 64

/**
 * \file mesh.h
 * \author Ryan Daly
 * \em rpdaly@ucsd.edu
 * \date 03/04/2009
 * This file contains the definition of the mesh_t struct which contains all mesh
 * information. It also contains constructors for all of the functions required to
 * manipulate meshes, including adding nodes and simplices.
 */

typedef struct
{
	int		    dim;			///< The Cartesian dimension in which the nodes lie
	long int	  nvert;			///< Number of vertices
	double		 *vertx;			///< List of nodes
	int		 *vtype;			///< List of node types
	long int	  nsimp;			///< Number of simplices
	long int	 *simpx;			///< Set of simplices
	int		 *ftype;			///< The face types; listed by combinadic index of vertex combinations from simplex
	long int	  nunff;			///< Number of unfinished faces
	long int	*ufaces;			///< List of unifinished faces
	long int	  nhulf;			///< Number of faces composing convex hull
	long int	*hfaces;			///< The faces belonging to the convex hull
	long int	**nsmap;			///< Node -> simplex map
	long int	*nssize;			///< nssize[i] = number of simplices to which vertex i belongs
} mesh_t;

/* Creating and destroying */
void	mesh_init( mesh_t *, long int );
void	mesh_free( mesh_t * );
int	mesh_load( mesh_t *, char *, int );		///< Loads a mesh from file; format specified by int argument
int	mesh_load_nodes( mesh_t *, char *, int );	///< Format specified by int
int	mesh_write( mesh_t *, char *, int );		///< Ditto
int	mesh_write_nodes( mesh_t *, char *, int );	///< Ditto

/* Important functions */
int	mesh_gen_node_simp_map( mesh_t * );		///< Fills in data needed by 
int	mesh_gen_init_simp( mesh_t *, long int );	///< This is a shitty function; don't use it unless absolutely necessary
int	mesh_gen_convex_hull( mesh_t * );		///< Generate the convex hull of faces on a vertex set
int	mesh_gen_del_wrap( mesh_t * );			///< Create a Delaunay mesh from a set of vertices
double	mesh_quality( mesh_t * );			///< Calculate the quality of the mesh
int	mesh_check_topology( mesh_t * );		///< Check for bad elements

/* Temporary output methods for debugging */
void	mesh_print_nodes( mesh_t * );
void	mesh_print_ufaces( mesh_t * );
void	mesh_print_nsmap( mesh_t * );

/* FIXME: TRY TO GET RID OF THE FUNCTIONS BELOW: THEY ARE NOT NEEDED. REPLACE THEM WITH list_t FUNCTIONS */
int	 mesh_del_unf_face( mesh_t *, long int );	///< Delete a face from ufaces at given position
int	 mesh_app_unf_face( mesh_t *, long int * );	///< Appends a face to ufaces
long int mesh_ufaces_index( mesh_t *, long int * );	///< Returns index of face of interest
int	 mesh_del_simpx( mesh_t *, long int );
int	 mesh_app_simpx( mesh_t *, long int * );
long int mesh_simpx_index( mesh_t *, long int * );
int	 mesh_del_hull_face( mesh_t *, long int );
int	 mesh_app_hull_face( mesh_t *, long int * );
long int mesh_hfaces_index( mesh_t *, long int * );
int	 mesh_set_hface_type( mesh_t *, long int, int );
int	 mesh_set_vertx_type( mesh_t *, long int, int );

/* FIXME: TRY TO GET RIGHT OF FUNCTIONS RETURNING int */
int	 mesh_vert_surf_test( mesh_t *, double *, long int * );			/* FIXME: GET RID OF ME		*/
int	 mesh_vert_surf_test_dim( long int, double *, double * );		/* FIXME: GET RID OF ME		*/
double	 mesh_vert_surf_test_interp_face( mesh_t *, double *, long int );	/* FIXME: RENAME THIS FUNCTION	*/
double	 mesh_vert_surf_test_interp_vert( mesh_t *, double *, long int * );	/* FIXME: RENAME THIS ONE TOO	*/
double	 mesh_face_parity( long int, double *, double * );
double	 mesh_simplex_determ( long int, double *, double * );
double	 mesh_vert_in_circsph_vert( mesh_t *, double *, long int * );		/* FIXME: RENAME THIS ONE TOO	*/
double	 mesh_vert_in_circsph( mesh_t *, double *, long int );			/* FIXME: RENAME THIS ONE TOO	*/

#endif
