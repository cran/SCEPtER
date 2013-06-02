#include <R.h>
#include <Rdefines.h>
#include<stdlib.h>
#include<Rmath.h>

typedef struct {
  double L, M, R;
} DATA;



SEXP lik2(SEXP data, SEXP star, SEXP sigma, SEXP thr, SEXP var)
{
  double *Pdata, *Psigma, *Pstar, *Pres;
  double *Teff, *logg, *z, *M, *R, *Dni, *nimax;
  double Vthr, maxL;
  long nrow, ncol, count;

  double sq2pi, chi[5], locsigma[5], chi2, mult, L, mass, radius;
  SEXP res, sel;
  long j, nres, DIM, n;
  int ii, norun, *Psel, *Pvar;
  DATA *d;

  // cast of input arguments and pointers
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
  sq2pi = sqrt(atan(1)*8);

  // dimensions of the dataset
  nrow = INTEGER(GET_DIM(data))[0];
  ncol = INTEGER(GET_DIM(data))[1];

  // column pointes
  Teff = Pdata;
  logg = Pdata+nrow;
  z = Pdata+2*nrow;
  Dni = Pdata+3*nrow;
  nimax = Pdata+4*nrow;
  M = Pdata+5*nrow;
  R = Pdata+6*nrow;

  // index vector for likelihood computations
  // 1 = include; 0 = exclude
  PROTECT(sel = NEW_INTEGER(nrow));
  Psel = INTEGER_POINTER(sel);
  for(j=0;j<nrow;j++)
    Psel[j] = 0;

  // compute mult for likelihood
  // scale sigma for Dni e nimax (in input is %)
  for(n=0;n<5;n++)
    locsigma[n] = Psigma[n];
  for(n=3;n<5;n++)
    locsigma[n] *= Pstar[n];
  
  mult = 1;
  for(n=0;n<5;n++)
    if(Pvar[n] == 1)
      mult *= 1.0/(sq2pi * locsigma[n]);

  // scan data and compute sel
  nres = 0;
  for(j=0;j<nrow;j++)
    {
      // chi
      for(ii=0;ii<5;ii++)
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

      norun = 0;
      for(ii=0;ii<5;ii++)
	{
	  if(fabs(chi[ii]) >= Vthr)
	    {
	      norun = 1;
	      break;
	    }
	}
      
      if( norun == 0 ) 
	{
	  nres++;
	  Psel[j] = 1;
	}
    }
  
  // no values! Return NULL
  if(nres == 0) 
    {
      UNPROTECT(6);
      return(R_NilValue);
    }

  // init the output matrix
  DIM = nres;
  d = (DATA *)calloc(DIM+1, sizeof(DATA));

  // compute lik, only for sel = 1
  nres = 0;
  maxL = 0;
  for(j=0;j<nrow;j++)
    {
      if( Psel[j] == 1 ) 
	{
	  // chi
	  for(ii=0;ii<5;ii++)
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

	  chi2 = 0;
	  for(n=0;n<5;n++)
	    chi2 += chi[n]*chi[n];
	  
	  L = mult * exp(-0.5*chi2);
	  if(L > maxL)
	    maxL = L;
	  d[nres].L = L;
	  d[nres].M = M[j];
	  d[nres].R = R[j];
	  nres++;
	}
    }

  mass = radius = 0;
  count = 0;
  for(j=0;j<nres;j++)
    {
      // select cases with lik >= 0.95 Max Lik
      if(d[j].L >= 0.95*maxL) 
	{
	  mass += d[j].M;
	  radius += d[j].R;
	  count++;
	}
    }
  mass /= (double)(count);
  radius /= (double)(count);
  
  // fill the output 
  PROTECT( res = NEW_NUMERIC(3) );
  Pres = NUMERIC_POINTER(res);
  Pres[0] = mass;
  Pres[1] = radius;
  Pres[2] = (double)count;

  free(d);
  
  // exit
  UNPROTECT(7);
  return(res);
}
