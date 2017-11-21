/* coalescentMCMC.c (2013-08-13) */

/* Copyright 2012-2013 Emmanuel Paradis

/* This file is part of the R-package `coalescentMCMC'. */
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

static R_CMethodDef C_entries[] = {
    {"get_single_index_integer", (DL_FUNC) &get_single_index_integer, 3},
    {"get_two_index_integer", (DL_FUNC) &get_two_index_integer, 3},
    {NULL, NULL, 0}
};

void R_init_coalescentMCMC(DllInfo *info)
{
    R_registerRoutines(info, C_entries, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
