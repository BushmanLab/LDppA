#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>

/* utils */

static void matmult(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;

  F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
		  x, &nrx, y, &nry, &zero, z, &nrx);
}

// //* As per Equation \ref{eq:prz}
void probzv( double *V, double *Z, int *T){
  double usable=1.0;
  for (int t = 0; t < *T; t++) {
    Z[t]=V[t]*usable;
    usable*=1-V[t];
  }
}
void zysum(
	   double *prw,
	   double *pz,
	   double *eps,
	   double *eta,
	   double *omdp,
	   int *w,
	   int *wp,
	   int *n,
	   int *T,
	   int *ka,
	   int *ko,
	   int *ndat,
	   /* result */
	   int *zy,
	   /* workspace */
	   double *etaomdp,
	   double *workT
	   ){
#ifdef TESTING 
  GetRNGstate();
#endif

  for (int idat=0;idat<*ndat;idat++){
    double biglog=R_NegInf;
    for (int t=0;t<*T;t++){
      int one=1L;
      // //* As per Equation \ref{eq:prwz3}. 
      // Note: eps==0 is limiting case 
      double x = (eps[0]==0) ? 0.0 :
	dnbinom(wp[t],eps[0],eps[1]/(eps[1]+prw[t]),one);
      for (int k=0;k<*ko;k++) x+= log(etaomdp[t+k**T]/prw[t]) *
				w[idat+k**ndat];
      if (biglog<x) biglog=x;
      workT[t]=x;
    }
    double prTot=0.0;
    for (int t=0;t<*T;t++){
      workT[t]=exp(workT[t]-biglog)*pz[t];
      prTot+=workT[t];
    }
    // //* As per multinomial part of \ref{eq:prw.z}
    for (int t=0;t<*T;t++) workT[t]/=prTot;
    rmultinom((int) n[idat],workT,(int) T[0], zy+idat**T);
  }  
#ifdef TESTING
  PutRNGstate();
#endif
}
/* main entry point                    */
/* run reps iterations of the sampler */

void innerLoop(
	       int *reps,         // nthin iterations
	       int *T,            // number of compositions
	       int *ka,           // number of actual cell types
	       int *ko,           // number of observed cell types
	       int *n,            // repetitions of a data pattern
	       int *ndat,         // number of data patterns
	       int *w,            // ndat by ko table of counts
	       int *wp,           // ndat vector of number of cells
	       double *s,         // prior for alpha
	       double *lamb,      // prior for row of eta
	       double *omcp,      // omega times 1-psi
	       double *omdp,      // omega times diag(psi)
	       double *eps,       // prior for R
	       /* modifiable */  
	       double *V,         // stick lengths
	       double *alpha,     // prior for V
	       double *eta,       // compositions
	       int *zy,           // T by ndat table of assignments
	       double *eoy,       // workspace
	       double *etaomdp,   // eta times omdp
	       double *prw,       // prob of observed count
	       double *pz,        // prob Z = t
	       double *workT,     // workspace
	       int *xstmp,        // workspace 
	       int *xsums         // T by ka table of counts
	       ){

  /* inits */
  GetRNGstate();

  for (int irep=0;irep<*reps;irep++){
    /* calc etaomdp and prw */
    matmult(eta, *T, *ka, omdp, *ka, *ko, etaomdp);

    /* rowsums */
    for (int t=0;t<*T;t++){
      double x=0.0;
      for (int k=0;k<*ka;k++) x+=etaomdp[t+k**T];
      prw[t]=x;
    }

    /* probzv */

    probzv(V, pz, T);

    /* zy */

    zysum(prw, pz, eps, eta, omdp, w, wp, n, T, ka, ko, ndat, zy, etaomdp, workT);

    /* xsums */

    for (int t=0;t<*T;t++){
      double seen=0.0;
      // //* See discussion of \ref{eq:xobs} 
      for (int k=0;k<*ka;k++) xsums[t+k**T]=0;
      for (int k=0;k<*ko;k++) {
        /* yz is number of cells in cluster t, cell type k */
        double yz=0.0;
        for (int idat=0;idat<*ndat;idat++){
          yz+=(double) (zy[t+idat**T]*w[idat+k**ndat]);
        }
        /* eoy: conditional prob of actual given observed */
        double eoysum=0.0;
        for (int k2=0;k2<*ka;k2++) {
          eoy[k2]=eta[t+k2**T]*omdp[k2+k**ka];
          eoysum+=eoy[k2];
        }
        for (int k2=0;k2<*ka;k2++) eoy[k2]/=eoysum;
        /* sample actual given observed */
        rmultinom( (int) yz, eoy, *ka, xstmp);
        for (int k2=0;k2<*ka;k2++) xsums[t+k2**T]+=xstmp[k2];
      }
      // //* Starting at \ref{eq:unseen}, see discussion
      for (int idat=0;idat<*ndat;idat++)
        seen+= (double) zy[t+idat**T] * (eps[0] + (double) wp[idat]);
      double nbparm = fmax2( 0.0, (1.0-prw[t])/(1.0+eps[1]) );
      /* R uses 1-nbparm in rnbinom where others might use nbparm */
      int notseen = (seen==0.0)? 0L: (int) rnbinom(seen,1.0-nbparm);
      /* Rprintf("seen = %f notseen = %d\n", seen, notseen); */
      if (notseen!=0){
        /* conditional prob of actual given unseen */
        double eoysum2=0.0;
        for (int k=0;k<*ka;k++){
          eoy[k]=eta[t+k**T]*omcp[k];
          eoysum2+=eoy[k];
        }
        for (int k=0;k<*ka;k++) eoy[k]/=eoysum2;
        /* sample actual given number unseen */
        rmultinom(notseen, eoy, *ka, xstmp);
        for (int k=0;k<*ka;k++) xsums[t+k**T]+=xstmp[k];
      }
      

    }

    /* eta */

    /* sample eta */
    // //* As per \ref{eq:eta.xl}
    for (int t=0; t<*T; t++){
      double x=0.0;
      for (int k=0;k<*ka;k++){
	int indx=t+k**T;
	eta[indx]=rgamma(xsums[indx]+lamb[k], 1.0);
	x+=eta[indx];
      }
      for (int k=0;k<*ka;k++){
	int indx=t+k**T;
	eta[indx]=eta[indx]/x;
      }
    }
    
    /* V */

    // //* See Equations \ref{eq:V.Za}
    V[*T-1L]=1.0;
    for (int t=0;t<*T;t++){
      int x=0.0;
      for (int idat=0;idat<*ndat;idat++)
        x+=zy[t+idat**T];
      workT[t]= (double) x;
    }
    double zgt=*alpha+workT[*T-1];
    for (int t=*T-2;t>=0;t--){
      double beta=rbeta(1.0+workT[t],zgt);
      V[t]= beta>0.99 ? 0.99 : beta;
      zgt+=workT[t];
    }

    /* alpha */

    // //* As per \ref{eq:alpha.r}
    double sumclog=0.0;
    for (int t=0; t<*T-1; t++) sumclog+=log(1.0-V[t]);
    alpha[0]=rgamma(s[0]+*T-1,1.0/(s[1]-sumclog));
  }
  /* fini */

  PutRNGstate();

}
