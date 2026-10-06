#ifndef SEMSOLVE2_SEMERROR_H
#define SEMSOLVE2_SEMERROR_H

#include <stdio.h>
#include <string.h>
#include <strings.h>

#define SEMERROR_MSGLEN 128

typedef struct
{
	int	id;
	int	line;
	char	file[SEMERROR_MSGLEN];
	char	message[SEMERROR_MSGLEN];
	char	function[SEMERROR_MSGLEN];
} semerror_t;

void semerror_full( semerror_t * );
void semerror_no_id( semerror_t * );
void semerror_no_message( semerror_t * );
void semerror_just_id( semerror_t * );

#endif
