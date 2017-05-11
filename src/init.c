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

extern void sampleCTC(void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *);

extern void zmat( void *prw, void *pz, void *eps, void *V, void *eta,
		  void *omdp, void *w, void *wp, void *n, void *T,
		  void *ka, void *ko, void *ndat, void *zy,
		  void *etaomdp, void *workT);

static const R_CMethodDef CEntries[] = {
  {"fixedEta",  (DL_FUNC) &fixedEta,  24},
  {"innerLoop", (DL_FUNC) &innerLoop, 24},
  {"sampleCTC", (DL_FUNC) &sampleCTC, 25},
  {"zmat",     (DL_FUNC)  &zmat,      16},
  {NULL, NULL, 0}
};

void R_init_LDppA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
