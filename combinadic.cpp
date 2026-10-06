#include <combinadic.h>

#include <cstdlib>
#include <cstring>
#include <strings.h>

/*
 * GENERAL COMBINATORIAL FUNCTIONS
 */

int factorial( int num_in )
{
	int fct = 1;
	int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

long int factorial( long int num_in )
{
	long int fct = 1;
	long int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

int binomial_slow( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return factorial( n_in ) / factorial( k_in ) / factorial( n_in - k_in );
}

long int binomial_slow( long int n_in, long int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return factorial( n_in ) / factorial( k_in ) / factorial( n_in - k_in );
}

int binomial( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact( n_in, k_in ) / factorial( k_in );
}

long int binomial( long int n_in, long int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact( n_in, k_in ) / factorial( k_in );
}

/**
 * This function returns the number of r_in-permutations of n_in objects.
 */
int truncfact( int n_in, int r_in )
{
	int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

long int truncfact( long int n_in, long int r_in )
{
	long int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

/**
 * Given a combinadic index, its corresponding vector is
 * the strictly increasing non-negative integer sequence,
 * {m_k}, of length dim_in and whose entries are less than
 * or equal to lim_in, and such that
 *   C(m_k,k+1) + C(m_{k-1},k) + ... + C(m_0,1) = idx_in
 */
int combinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	int i,j,idx,tmp;

	if( idx_in < 0 )
		return -1;

	idx = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		for(j=i+1;j<=lim_in;j++)
		{
			tmp = binomial( j, i + 1 );
			if( tmp > idx )
				break;
		}
		idx -= binomial( j - 1, i + 1 );
		vec_out[i] = j - 1;
	}
	return 0;
}

long int combinadic_vector( long int idx_in, long int lim_in, long int dim_in, long int *vec_out )
{
	long int i,j,idx,tmp;

	if( idx_in < 0 )
		return -1;

	idx = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		for(j=i+1;j<=lim_in;j++)
		{
			tmp = binomial( j, i + 1 );
			if( tmp > idx )
				break;
		}
		idx -= binomial( j - 1, i + 1 );
		vec_out[i] = j - 1;
	}
	return 0;
}

int combinadic_vector_length( int idx_in, int dim_in, int lim_in )
{
	int i,sum;

	sum = 0;
	for(i=1;i<=dim_in;i++)
	{
		sum += binomial( lim_in, i );
		if( sum > idx_in )
			break;
	}
	return i;
}

long int combinadic_vector_length( long int idx_in, long int dim_in, long int lim_in )
{
	long int i,sum;

	sum = 0;
	for(i=1;i<=dim_in;i++)
	{
		sum += binomial( lim_in, i );
		if( sum > idx_in )
			break;
	}
	return i;
}

int combinadic_next( int lim_in, int dim_in, int *vec_in )
{
	int i,j;
	for(i=0;i<dim_in;i++)
	{
		if( ( vec_in[i] < lim_in )
			&& ( ( i < dim_in - 1 && ( vec_in[i+1] - vec_in[i] ) > 1 ) || i == dim_in - 1 ) )
		{
			++vec_in[i];
			break;
		}
	}
	for(j=0;j<i;j++)
		vec_in[j] = j;
}

long int combinadic_next( long int lim_in, long int dim_in, long int *vec_in )
{
	long int i,j;
	for(i=0;i<dim_in;i++)
	{
		if( ( vec_in[i] < lim_in )
			&& ( ( i < dim_in - 1 && ( vec_in[i+1] - vec_in[i] ) > 1 ) || i == dim_in - 1 ) )
		{
			++vec_in[i];
			break;
		}
	}
	for(j=0;j<i;j++)
		vec_in[j] = j;
}

int combinadic_init( int lim_in, int dim_in, int *ptr_in )
{
	int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

long int combinadic_init( long int lim_in, long int dim_in, long int *ptr_in )
{
	long int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

int combinadic_index( int dim_in, int *vec_in )
{
	int i;
	int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}

long int combinadic_index( long int dim_in, long int *vec_in )
{
	long int i;
	long int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}

#ifdef __cplusplus
int combinadic_index( int dim_in, std::vector< int > &vec_in )
{
	int i;
	int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}

long int combinadic_index( long int dim_in, std::vector< long int > &vec_in )
{
	long int i;
	long int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}
#endif

/*
 * FACTORADIC OPERATIONS
 */

int factoradic_radix_index( int dim_in, int idx_in, int *rad_out )
{
	int i,fct,div,rem;

	rem = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		fct = factorial( i );
		div = rem / fct;
		rem = rem % fct;
		rad_out[dim_in-1-i] = div;
	}
}

int factoradic_radix_vector( int dim_in, int *vec_in, int *rad_out )
{
	int c,i,j;
	for(i=0;i<dim_in;i++)
		rad_out[i] = 0;
	for(i=0;i<dim_in;i++)
	{
		for(j=i+1;j<dim_in;j++)
		{
			if( vec_in[j] < vec_in[i] )
				++rad_out[i];
		}
	}
	return 0l;
}

int factoradic_index( int dim_in, int *vec_in )
{
	int i,rad[dim_in],idx = 0;
	factoradic_radix_vector( dim_in, vec_in, rad );
	for(i=0;i<dim_in;i++)
		idx += rad[dim_in-1-i] * factorial( i );
	return idx;
}

int factoradic_init( int dim_in, int *vec_in )
{
	int i;

	for(i=0;i<dim_in;i++)
		vec_in[i] = i;

	return 0l;
}

int factoradic_next( int dim_in, int *vec_in )
{
	
}

int factoradic_vector( int idx_in, int dim_in, int *vec_out )
{
	int i,j,k,rad[dim_in];
	int fnd[dim_in];

	for(i=0;i<dim_in;i++)
		fnd[i] = 0;

	factoradic_radix_index( dim_in, idx_in, rad );

	for(i=0;i<dim_in;i++)
	{
		for(j=0,k=0;j<dim_in;j++)
		{
			if( fnd[j] == 0 )
				++k;
			if( k - 1 >= rad[i] )
				break;
		}
		fnd[j] = 1;
		vec_out[i] = j;
	}
	return 0;
}

/**
 * Given a factoradic index, this function computes the radix form of the
 * index. That is, it converts the base 10 input index into a non-constant
 * factorial base number.
 */
long int factoradic_radix_index_long( long int dim_in, long int idx_in, long int *rad_out )
{
	long int i,fct,div,rem;

	rem = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		fct = factorial( i );
		div = rem / fct;
		rem = rem % fct;
		rad_out[dim_in-1-i] = div;
	}
}

long int factoradic_radix_vector_long( long int dim_in, long int *vec_in, long int *rad_out )
{
	long int c,i,j;
	for(i=0;i<dim_in;i++)
		rad_out[i] = 0;
	for(i=0;i<dim_in;i++)
	{
		for(j=i+1;j<dim_in;j++)
		{
			if( vec_in[j] < vec_in[i] )
				++rad_out[i];
		}
	}
	return 0l;
}

long int factoradic_index_long( long int dim_in, long int *vec_in )
{
	long int i,rad[dim_in],idx = 0;
	factoradic_radix_vector_long( dim_in, vec_in, rad );
	for(i=0;i<dim_in;i++)
		idx += rad[dim_in-1-i] * factorial( i );
	return idx;
}

long int factoradic_init_long( long int dim_in, long int *vec_in )
{
	long int i;

	for(i=0;i<dim_in;i++)
		vec_in[i] = i;

	return 0l;
}

long int factoradic_next_long( long int dim_in, long int *vec_in )
{
	
}

long int factoradic_vector_long( long int idx_in, long int dim_in, long int *vec_out )
{
	long int i,j,k,rad[dim_in];
	int fnd[dim_in];

	for(i=0;i<dim_in;i++)
		fnd[i] = 0;

	factoradic_radix_index_long( dim_in, idx_in, rad );

	for(i=0;i<dim_in;i++)
	{
		for(j=0,k=0;j<dim_in;j++)
		{
			if( fnd[j] == 0 )
				++k;
			if( k - 1 >= rad[i] )
				break;
		}
		fnd[j] = 1;
		vec_out[i] = j;
	}
	return 0;
}

/*
 * BINARY SEARCH ROUTINES
 */

int binary_search( long int *i_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	left = 0;
	right = *llen_in - 1;
	return binary_search_range( i_in, &left, &right, list_in, llen_in, idx_out, its_out );
}

int binary_search_range( long int *i_in, long int *left_in, long int *right_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	*its_out = 0;
	if( *llen_in > 0 )
	{
		left = *left_in;
		right = *right_in;
		while(left<=right)
		{
			*idx_out = ( left + right ) / 2;
			if( list_in[*idx_out] == *i_in )
				return 0;
			else
			{
				if( list_in[*idx_out] < *i_in )
					left = *idx_out + 1;
				else
					right = *idx_out - 1;
			}
			*its_out = *its_out + 1;
		}
		return -1;
	}
	else
	{
		*idx_out = 0;
		return -1;
	}
}

int binary_closest( long int *i_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	left = 0;
	right = *llen_in - 1;
	return binary_closest_range( i_in, &left, &right, list_in, llen_in, idx_out, its_out );
}

int binary_closest_range( long int *i_in, long int *left_in, long int *right_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int idx;
	if( binary_search_range( i_in, left_in, right_in, list_in, llen_in, &idx, its_out ) == -1 )
	{
		if( *llen_in > 0 && idx > -1 )
			*idx_out = list_in[idx] > *i_in ? idx - 1 : idx;
		else
			*idx_out = -1;
		return -1;
	}
	else
	{
		*idx_out = idx;
		return 0;
	}
}

/**
 * Just a simple search-for-integer function. Returns -1 if not found.
 */
int simple_find( int target_in, int len_in, int *set_in )
{
	int i;
	for(i=0;i<len_in;i++)
	{
		if( set_in[i] == target_in )
			return i;
	}
	return -1;
}

int simple_is_subset( int len1_in, int *set1_in, int len2_in, int *set2_in )
{
	int i,j;
	for(i=0;i<len1_in;i++)
	{
		if( simple_find( set1_in[i], len2_in, set2_in ) == -1 )
			return 0; // set1_in is NOT a subset of set2_in
	}
	return 1; // set1_in is a subset of set2_in
}

int simple_set_intersect( int len1_in, int *set1_in, int len2_in, int *set2_in, int *set_out )
{
	int i,j;
	for(i=0,j=0;i<len1_in;i++)
	{
		if( simple_find( set1_in[i], len2_in, set2_in ) != -1 )
		{
			set_out[j] = set1_in[i];
			++j;
		}
	}
	return j;
}

void bubble_sort( int len_in, int *set_in )
{
	int i,j,tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

/**
 * Just a simple search-for-integer function. Returns -1 if not found.
 * This one is just the long version.
 */
long int simple_find_long( long int target_in, long int len_in, long int *set_in )
{
	long int i;
	for(i=0;i<len_in;i++)
	{
		if( set_in[i] == target_in )
			return i;
	}
	return -1;
}

long int simple_is_subset_long( long int len1_in, long int *set1_in, long int len2_in, long int *set2_in )
{
	long int i;
	for(i=0;i<len1_in;i++)
	{
		if( simple_find_long( set1_in[i], len2_in, set2_in ) == -1 )
			return 0; // set1_in is NOT a subset of set2_in
	}
	return 1; // set1_in is a subset of set2_in
}

long int simple_set_intersect_long( long int len1_in, long int *set1_in, long int len2_in, long int *set2_in, long int *set_out )
{
	long int i,j;
	for(i=0,j=0;i<len1_in;i++)
	{
		if( simple_find_long( set1_in[i], len2_in, set2_in ) != -1 )
		{
			set_out[j] = set1_in[i];
			++j;
		}
	}
	return j;
}

void bubble_sort_long( long int len_in, long int *set_in )
{
	long int i,j,tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

void bubble_sort_double_long( long int len_in, double *set_in )
{
	long int i,j;
	double tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

int list_init( list_t *list_in, int step_in )
{
	void *tmp;

	list_in->start = NULL;
	list_in->pointer = NULL;
	list_in->step = step_in;
	list_in->size = 0;
	list_in->alloc = LIST_SIZE_INC;
	tmp = std::realloc( list_in->pointer, LIST_SIZE_INC * step_in );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	return 0;
}

int list_free( list_t *list_in )
{
	if( list_in->start != NULL && list_in->size > 0 )
		free( list_in->pointer );
	else
		return -1;
	list_in->pointer = NULL;
	list_in->start = NULL;
	return 0;
}

int list_reset( list_t *list_in )
{
	if( list_in->size <= 0 )
		return -1;
	std::memmove( list_in->pointer, list_in->start, list_in->size * list_in->step );
	list_in->start = list_in->pointer;
	return 0;
}

int list_shift( list_t *list_in, int shift_in )
{
	if( list_in->size <= 0
	 || ( ( (byte_t*) list_in->start + shift_in * list_in->step ) < list_in->pointer )
	 || ( ( (byte_t*) list_in->start + ( list_in->size + shift_in ) * list_in->step ) > ( (byte_t*) list_in->pointer + list_in->alloc * list_in->step ) ) )
		return -1;
	if( shift_in == 0 )
		return 0;
	list_in->start = std::memmove( (byte_t*) list_in->start + shift_in * list_in->step, list_in->start, list_in->size * list_in->step );
	return 0;
}

int list_grow( list_t *list_in, int minc_in )
{
	void *tmp;
	list_reset( list_in );
	tmp = realloc( list_in->pointer, ( list_in->alloc + minc_in ) * list_in->step );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc += minc_in;
	return list_in->alloc;
}

int list_shrink( list_t *list_in, int mdec_in )
{
	void *tmp;
	list_reset( list_in ); /* Force list_in->start to coincide with list_in->pointer */
	tmp = realloc( list_in->pointer, ( list_in->alloc - mdec_in ) * list_in->step );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc -= mdec_in;
	return list_in->alloc;
}

int list_prepend( list_t *list_in, void *val_in )
{
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
		if( list_shift( list_in, 1 ) == -1 )
			return -1;
	}
	else
	{
		if( list_in->start <= list_in->pointer ) /* Meaning there is no memory allocated to the left, so make some */
			if( list_shift( list_in, 1 ) == -1 )
				return -1;
	}
	list_in->start = (byte_t*) list_in->start - list_in->step;
	bcopy( val_in, list_in->start, list_in->step );
	++list_in->size;
	return list_in->size;
}

int list_append( list_t *list_in, void *val_in )
{
	/* If there isn't enough memory allocated */
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	}
	else /* If there is enough, then either there is some at the beginning or some at the end */
	{
		if( (byte_t*)list_in->start > (byte_t*)list_in->pointer )
			if( list_shift( list_in, -1 ) == -1 )
				return -1;
	}

	bcopy( val_in, (byte_t*) list_in->start + list_in->size * list_in->step, list_in->step );
	++list_in->size;
	return list_in->size;
}

int list_popf( list_t *list_in, void *ptr_out )
{
	if( list_in->size <= 0 )
		return -1;
	bcopy( list_in->start, ptr_out, list_in->step );
	if( list_in->size > 1 )
		list_in->start = (byte_t*) list_in->start + list_in->step;
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		if( list_shrink( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	return 0;
}

int list_popb( list_t *list_in, void *ptr_out )
{
	if( list_in->size <= 0 )
		return -1;
	else
		bcopy( (byte_t*) list_in->start + ( list_in->size - 1 ) * list_in->step, ptr_out, list_in->step );
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		if( list_shrink( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	return 0;
}
