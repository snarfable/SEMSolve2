%{
	#include <stdio.h>
	#include <stdlib.h>
	#include "parse_symb.h"
%}

%%
commands:
	| commands command
	;
command:
	
	;
%%

/* Define yyerror at the bottom of your .y file */
void yyerror(const char *s) {
    fprintf(stderr, "Error: %s\n", s);
}
