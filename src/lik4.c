#include <R.h>
#include <Rdefines.h>
#include<stdlib.h>
#include<Rmath.h>

#define NVAR 7  

typedef struct {
  double L, M, R, logage;
} DATA;


void findrange(double *x, long dim, double vl, double vu, long *l, long *u)
{
  long res, up, low;
  
  if(vl <= x[0])
    res = 0;
  else if(vl > x[dim-1])
    res = -1;
  else
    {
      up = dim-1;
      low = 0;
      while(up-low > 1)
	{
	  res = (up+low)/2;
	  if(x[res] < vl)
	    low = res;
	  else 
	    up = res;
	}
      res = up;
    }
  *l = res;

  if(vu < x[0])
    res = -1;
  else if(vu >= x[dim-1])
    res = dim-1;
  else
    {
      up = dim-1;
      low = *l;
      while(up-low > 1)
	{
	  res = (up+low)/2;
	  if(x[res] < vu)
	    low = res;
	  else 
	    up = res;
	}
      res = low;
    }
  *u = res;
}


SEXP lik4(SEXP data, SEXP star, SEXP sigma, SEXP thr, SEXP var)
{
  double *Pdata, *Psigma, *Pstar, *Pres;
  double *Teff, *logg, *z, *M, *R, *Dni, *nimax, *logage;
  double Vthr, maxL, lmult;
  long nrow, ncol, count;

  double sq2pi, chi[NVAR], locsigma[NVAR], chi2, mult, L, mass, 
    radius, lt, ltnlog;;
  double sTeffP, sTeffM;
  SEXP res, dm, sel;
  long i, j, nres, DIM, start, n, startT, stopT, up, low;
  int ii, norun, *Psel, *Pvar;
  DATA *d;

  // cast degli argomenti e creazione loro puntatori
  PROTECT(data = AS_NUMERIC(data));
  PROTECT(star = AS_NUMERIC(star));
  PROTECT(sigma = AS_NUMERIC(sigma));
  PROTECT(thr = AS_NUMERIC(thr));
  PROTECT(var = AS_INTEGER(var));

  Pdata = NUMERIC_POINTER(data);
  Pstar = NUMERIC_POINTER(star);
  Psigma = NUMERIC_POINTER(sigma);
  Vthr = NUMERIC_VALUE(thr);
  Pvar = INTEGER_POINTER(var);
  
  // sqrt ( 2 * pi )
  sq2pi = 2.506628274631000;

  // dataset dimensions
  nrow = INTEGER(GET_DIM(data))[0];
  ncol = INTEGER(GET_DIM(data))[1];

  // column pointers
  Teff = Pdata;
  logg = Pdata+nrow;
  z = Pdata+2*nrow;
  Dni = Pdata+3*nrow;
  nimax = Pdata+4*nrow;
  M = Pdata+5*nrow;
  R = Pdata+6*nrow;
  logage = Pdata+7*nrow;

  // index vector for likelihood computations
  // 1 = include; 0 = exclude
  Psel = (int*)malloc(nrow*sizeof(int));
  for(j=0;j<nrow;j++)
    Psel[j] = 0;

  // compute mult for likelihood
  // scale sigma for Dni e nimax (in input is %)
  for(n=0;n<NVAR;n++)
    locsigma[n] = Psigma[n];
  for(n=3;n<5;n++)
    locsigma[n] *= Pstar[n];
  
  mult = 1;
  for(n=0;n<NVAR;n++)
    if(Pvar[n] == 1)
      mult *= 1.0/(sq2pi * locsigma[n]);
  lmult = log(mult);

  // Teff range for selection
  sTeffP = Pstar[0] + Vthr*locsigma[0];
  sTeffM = Pstar[0] - Vthr*locsigma[0];


  // the dataset must be ordered by Teff !!
  findrange(Teff, nrow, sTeffM, sTeffP, &startT, &stopT);
  if(startT == -1 || stopT == -1) {
    UNPROTECT(5);
    return(R_NilValue);
  }
  
  // scan data and compute sel
  nres = 0;
  for(j=startT;j<=stopT;j++)
    {
      // chi
      for(ii=0;ii<NVAR;ii++)
	chi[ii] = 0;

      if(Pvar[0] == 1)
	chi[0] = (Teff[j] - Pstar[0])/locsigma[0];
      if(Pvar[1] == 1)
	chi[1] = (logg[j] - Pstar[1])/locsigma[1];
      if(Pvar[2] == 1)
	chi[2] = (z[j] - Pstar[2])/locsigma[2];
      if(Pvar[3] == 1)
	chi[3] = (Dni[j] - Pstar[3])/locsigma[3];
      if(Pvar[4] == 1)
	chi[4] = (nimax[j] - Pstar[4])/locsigma[4];
      if(Pvar[5] == 1)
	chi[5] = (M[j] - Pstar[5])/locsigma[5];
      if(Pvar[6] == 1)
	chi[6] = (R[j] - Pstar[6])/locsigma[6];

      norun = 0;
      for(ii=0;ii<NVAR;ii++)
	{
	  if(fabs(chi[ii]) >= Vthr)
	    {
	      norun = 1;
	      break;
	    }
	}
      
      if( norun == 0 ) 
	{
	  chi2 = 0;
	  for(ii=0;ii<NVAR;ii++)
	    chi2 += chi[ii]*chi[ii];
	  
	  nres++;
	  Psel[j] = 1;
	}
    }

  // no values! Return NULL
  if(nres == 0) 
    {
      free(Psel);
      UNPROTECT(5);
      return(R_NilValue);
    }

  // init the output matrix
  DIM = nres;
  d = (DATA *)calloc(DIM+1, sizeof(DATA));
  PROTECT( res = NEW_NUMERIC(5) );
  Pres = NUMERIC_POINTER(res);

  // compute lik, only for sel = 1
  nres = 0;
  maxL = -1e99;
  for(j=startT;j<=stopT;j++)
    {
      if( Psel[j] == 1 ) 
	{
	  // chi
	  for(ii=0;ii<NVAR;ii++)
	    chi[ii] = 0;
	  
	  if(Pvar[0] == 1)
	    chi[0] = (Teff[j] - Pstar[0])/locsigma[0];
	  if(Pvar[1] == 1)
	    chi[1] = (logg[j] - Pstar[1])/locsigma[1];
	  if(Pvar[2] == 1)
	    chi[2] = (z[j] - Pstar[2])/locsigma[2];
	  if(Pvar[3] == 1)
	    chi[3] = (Dni[j] - Pstar[3])/locsigma[3];
	  if(Pvar[4] == 1)
	    chi[4] = (nimax[j] - Pstar[4])/locsigma[4];
	  if(Pvar[5] == 1)
	    chi[5] = (M[j] - Pstar[5])/locsigma[5];
	  if(Pvar[6] == 1)
	    chi[6] = (R[j] - Pstar[6])/locsigma[6];

	  chi2 = 0;
	  for(n=0;n<NVAR;n++)
	    chi2 += chi[n]*chi[n];
	  
	  // log likelihood
	  L = lmult -0.5*chi2;
	  if(L > maxL)
	    maxL = L;
	  d[nres].L = L;
	  d[nres].M = M[j];
	  d[nres].R = R[j];
	  d[nres].logage = logage[j];
	  nres++;
	}
    }

  mass = radius = lt = ltnlog = 0;
  count = 0;
  maxL = log(0.95) + maxL;
  for(j=0;j<nres;j++)
    {
       // select cases with lik >= 0.95 Max Lik
      if(d[j].L >= maxL) 
	{
	  mass += d[j].M;
	  radius += d[j].R;
	  lt += d[j].logage;
	  ltnlog += 1e-9*pow(10, d[j].logage);
	  count++;
	}
    }
  mass /= (double)(count);
  radius /= (double)(count);
  lt /= (double)(count);
  ltnlog /= (double)(count);

  // fill the output 
  Pres[0] = mass;
  Pres[1] = radius;
  Pres[2] = lt;
  Pres[3] = ltnlog;
  Pres[4] = (double)count;

  free(d);
  free(Psel);
  
  // exit
  UNPROTECT(6);
  return(res);
}

