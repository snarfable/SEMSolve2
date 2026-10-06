%{
	#include <stdio.h>
	#include <stdlib.h>
	#include "parse_symb.h"

%}

%%
=			return EQUAL_T;
"\["			return LBRACKET_T;
"\]"			return RBRACKET_T;
[ \n\t]			/* Ignore white space */
#[^\n]+			/* Comments */
%%
