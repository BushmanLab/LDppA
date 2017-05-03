#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */

extern void fixedEta(void *, void *, void *, void *, void *,
		     void *, void *, void *, void *, void *,
		     void *, void *, void *, void *, void *,
		     void *, void *, void *, void *, void *,
		     void *, void *, void *, void *);

extern void innerLoop(void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"fixedEta",  (DL_FUNC) &fixedEta,  24},
  {"innerLoop", (DL_FUNC) &innerLoop, 24},
  {NULL, NULL, 0}
};

void R_init_LDppA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}