#ifndef COMBINADIC_H
#define COMBINADIC_H

#ifdef __cplusplus
#include <vector>
#endif

#define LIST_SIZE_INC 20
#define LIST_NULL_PTR 0xffffffffffffffff

/* Typical versions */
int      factorial( int );
long int factorial( long int );
int      binomial_slow( int n_in, int k_in );
long int binomial_slow( long int n_in, long int k_in );
int      binomial( int, int );
long int binomial( long int, long int );
int      truncfact( int, int );
long int truncfact( long int, long int );

/* Combinadics */
int combinadic_index( int, int * ); 				/* dim_in, vec_in, returns index	*/
#ifdef __cplusplus
int combinadic_index( int, std::vector< int > & );		/* dim_in, vec_in, returns index	*/
#endif
int combinadic_init( int, int, int * );				/* lim_in, dim_in, vec_out		*/
int combinadic_next( int, int, int * );				/* lim_in, dim_in, cur_vec_in/out	*/
int combinadic_vector( int, int, int, int * );			/* index_in, lim_in, dim_in, vec_out	*/
int combinadic_vector_length( int, int, int );

/* Long combinadics */
long int combinadic_index( long int, long int * );
#ifdef __cplusplus
long int combinadic_index( long int, std::vector< long int > & );
#endif
long int combinadic_init( long int, long int, long int * );
long int combinadic_next( long int, long int, long int * );
long int combinadic_vector( long int, long int, long int, long int * );
long int combinadic_vector_length( long int, long int, long int );

/* Factoradics */
int factoradic_radix_index( int, int, int * );
int factoradic_radix_vector( int, int *, int * );
int factoradic_index( int, int * );
int factoradic_init( int, int * );
int factoradic_next( int, int * );
int factoradic_vector( int, int, int * );

/* Long factoradics */
long int factoradic_radix_index_long( long int, long int, long int * );
long int factoradic_radix_vector_long( long int, long int *, long int * );
long int factoradic_index_long( long int, long int * );
long int factoradic_init_long( long int, long int * );
long int factoradic_next_long( long int, long int * );
long int factoradic_vector_long( long int, long int, long int * );

/* Stirling numbers of the second kind */
int stirling_index( int, int * ); ///< dimension, input vector
int stirling_vector( int, int, int * ); ///< index, dimension, output vector

/* Long Stirling numbers */
long int stirling_index_long( long int, long int * ); ///< dimension, input vector
long int stirling_vector_long( long int, long int, long int * ); ///< index, dimension, output vector

/* Binary searches */
int binary_search_int( int *, int *, int *, int *, int * );
int binary_search_range_int( int *, int *, int *, int *, int *, int *, int * );
int binary_closest_int( int *, int *, int *, int *, int * );
int binary_closest_range_int( int *, int *, int *, int *, int *, int *, int * );
int binary_set_intersect( int, int *, int, int *, int * );

/* Long versions */
int binary_search( long int *, long int *, long int *, long int *, long int * );
int binary_search_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * );
int binary_closest( long int *, long int *, long int *, long int *, long int * );
int binary_closest_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * );
int binary_set_intersect_long( long int, long int *, long int, long int *, long int * );

/* Simple searches and sorts */
int simple_find( int, int, int * );
int simple_is_subset( int, int *, int, int * );
int simple_set_intersect( int, int *, int, int *, int * );
void bubble_sort( int, int * );
void bubble_sort_double( int, double * );
long int simple_find_long( long int, long int, long int * );
long int simple_is_subset_long( long int, long int *, long int, long int * );
long int simple_set_intersect_long( long int, long int *, long int, long int *, long int * );
void bubble_sort_long( long int, long int * );
void bubble_sort_double_long( long int, double * );

typedef struct
{
	void *start;
	void *pointer;
	int step;
	int size;
	int alloc;
} list_t;

typedef char byte_t;

int list_init( list_t *, int );
int list_free( list_t * );
int list_reset( list_t * );
int list_shift( list_t *, int );
int list_grow( list_t *, int );
int list_shrink( list_t *, int );
int list_prepend( list_t *, void * );
int list_append( list_t *, void * );
int list_popf( list_t *, void * );
int list_popb( list_t *, void * );

#endif
