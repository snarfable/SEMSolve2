#include <semerror.h>

void semerror_full( semerror_t *err_in )
{
	printf( "[Error %d]: %s(): in file `%s' @ line %d: %s\n", err_in->id, err_in->function, err_in->file, err_in->line, err_in->message );
}

void semerror_no_id( semerror_t *err_in )
{
	printf( "[Error]: in file `%s' @ line %d: %s\n", err_in->file, err_in->line, err_in->message );
}

void semerror_no_message( semerror_t *err_in )
{
	printf( "[Error %d]: in file `%s' @ line %d\n", err_in->id, err_in->file, err_in->line );
}

void semerror_just_id( semerror_t *err_in )
{
	printf( "[Error %d]: %s\n", err_in->id, err_in->message );
}
