#ifndef SEMSOLVE2_SEMCONFIG_H
#define SEMSOLVE2_SEMCONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

typedef struct
{
	char		file[128];
	xmlDocPtr	docptr;
} semconfig_t;

#endif
