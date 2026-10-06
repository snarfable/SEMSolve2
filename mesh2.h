#ifndef SEMSOLVE2_MESH2_H
#define SEMSOLVE2_MESH2_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <combinadic.h>

typedef struct
{
	long int dim;
	long int nvert;
	long int nsimp;
	list_t *vertx;
	list_t *simpx;
	list_t *uface;
	list_t *nsmap;
} mesh_t;

/* Management functions */
int mesh_init( mesh_t *, long int );
int mesh_free( mesh_t * );
int mesh_load( mesh_t *, char *, int );
int mesh_write( mesh_t *, char *, int );
int mesh_load_nodes( mesh_t *, char *, int );
int mesh_write_nodes( mesh_t *, char *, int );

/* Important functions */
int mesh_gen_node_simp_map( mesh_t * );
int mesh_gen_convex_hull( mesh_t * );
int mesh_gen_del_wrap( mesh_t * );

/* Output functions for debugging */
void mesh_print_nodes( mesh_t * );
void mesh_print_ufaces( mesh_t * );
void mesh_print_nsmap( mesh_t * );

/* Functions for face checking */
double mesh_face_test( long int, double *, double * );
double mesh_face_test_interp( long int, double *, double * );
double mesh_face_parity( long int, double *, double * );
double mesh_circumsphere( long int, double *, double * );

#endif
