#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

/* .C calls */

extern void sampleCTC(void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *,
		      void *);

extern void zmat( void *prw, void *pz, void *eps, void *V, void *eta,
		  void *omdp, void *w, void *wp, void *n, void *T,
		  void *ka, void *ko, void *ndat, void *zy,
		  void *etaomdp, void *workT);

static const R_CMethodDef CEntries[] = {
  {"sampleCTC", (DL_FUNC) &sampleCTC, 26},
  {"zmat",     (DL_FUNC)  &zmat,      16},
  {NULL, NULL, 0}
};


/* .Call calls */

extern SEXP amllk( SEXP phi, SEXP w, SEXP omegaPsi);
extern SEXP dldphi( SEXP phi, SEXP w, SEXP op);
extern SEXP fixedEtaEMCall(SEXP lkMatrix, SEXP probg,
			   SEXP ritab, SEXP alpha,
			   SEXP weights, SEXP nreps);

static const R_CallMethodDef CallEntries[] = {
  CALLDEF(amllk, 3),
  CALLDEF(dldphi, 3),
  CALLDEF(fixedEtaEMCall,6),
  {NULL, NULL, 0}
};

void R_init_ECTC(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
