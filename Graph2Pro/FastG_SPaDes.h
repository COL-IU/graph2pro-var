#ifndef __FASTGSPADEX_H_
#define __FASTGSPADEX_H_

#include "lib.h"

namespace FastGSPaDes
{
	int countedge(char *edgefile);
	int inputedgenames(char *edgeseqfile, char **alledgename);
	int locateedge(char **alledgename, char *edgename, int num_edge);
	int locaterevedge(char **alledgename, char *edgename, int num_edge);
	int readedge(char *edgeseqfile, EDGE *edge, VERTEX *vertex, int *edgeindex, int num_edge, int kmersize, char **alledgename);
	int strcmpr(const void * a, const void * b);
	int get_next_edge(char *tmpstr, int *start, int end, char **alledgename, int num_edge);
	uint256_t str2code(char *seq, int kmersize);
};

#endif
