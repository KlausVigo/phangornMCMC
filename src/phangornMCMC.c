/* phangornMCMC.c (2013-08-13) */

/* Copyright 2012-2013 Emmanuel Paradis

/* This file copied of the R-package `coalescentMCMC'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <R_ext/Rdynload.h>

void get_single_index_integer(int *x, int *val, int *index)
{
	int i = 0, v = *val;
	while (x[i] != v) i++;
	*index = i + 1;
}

void get_two_index_integer(int *x, int *val, int *index)
{
	int i1 = 0, i2, v = *val;
	while (x[i1] != v) i1++;
	i2 = i1 + 1;
	while (x[i2] != v) i2++;
	index[0] = i1 + 1;
	index[1] = i2 + 1;
}


