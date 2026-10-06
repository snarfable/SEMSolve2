#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ltgmap.h>
#include <combinadic.h>

int main( int argc, char **argv )
{
	list_t list,temp;
	list_init( &list, sizeof(long int) );
	int i;
	long int tmp;
	double dbl;
	for(i=0;i<30;i++)
	{
		list_append( &list, &i );
		printf( "list.alloc = %d\n", list.alloc );
	}
	for(i=0;i<20;i++)
	{
		list_prepend( &list, &i );
		printf( "now list.alloc = %d\n", list.alloc );
	}
	for(i=0;i<5;i++)
	{
		list_popf( &list, &tmp );
		printf( "popf = %d\n", tmp );
		printf( "NOW list.alloc = %d, list.size = %d\n", list.alloc, list.size );
	}
	for(i=0;i<10;i++)
	{
		list_popb( &list, &tmp );
		printf( "popb = %d\n", tmp );
		printf( "NOW list.alloc = %d, list.size = %d\n", list.alloc, list.size );
	}
	for(i=0;i<16;i++)
	{
		list_popf( &list, &tmp );
		printf( "popf = %d\n", tmp );
		printf( "NOW list.alloc = %d, list.size = %d\n", list.alloc, list.size );
	}
	for(i=0;i<list.size;i++)
		printf( "%d = %d\n", i, ((long int*)list.start)[i] );

	printf( "list.start = %x\nlist.pointer = %x\nlist.size = %d\nlist.alloc = %d\n", list.start, list.pointer, list.size, list.alloc );
	list_free( &list );

	list_init( &list, sizeof(list_t) );
	list_init( &temp, sizeof(double) );
	for(i=0;i<5;i++)
		dbl = 2.58 / (double) i, list_append( &temp, &dbl );
	for(i=0;i<10;i++)
		list_append( &list, &temp );
	list_free( &temp );
	list_free( &list );

	return 0;
}
