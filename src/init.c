#include <R.h>
#include <R_ext/Rdynload.h>


/* .C calls */
extern void get_single_index_integer(void *, void *, void *);
extern void get_two_index_integer(void *, void *, void *);


static R_CMethodDef C_entries[] = {
    {"get_single_index_integer", (DL_FUNC) &get_single_index_integer, 3},
    {"get_two_index_integer",    (DL_FUNC) &get_two_index_integer,    3},
    {NULL, NULL, 0}
};

void R_init_phangornMCMC(DllInfo *info)
{
    R_registerRoutines(info, C_entries, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
