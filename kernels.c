#include <stdio.h>
#include <math.h>

#include "header.h"

float cubic_spline_kernel(float u)
{
  float prefac;
  
#if NDIM==3
  prefac=1./M_PI;
#elif NDIM==2
  prefac=10./7./M_PI;
#else
  prefac=2./3.;
#endif  

  if(u<=1)
    return prefac*(1.-1.5*u*u*(1.-0.5*u));
  else if(u>1 && u<2)
    return prefac*0.25*pow(2.-u,3.0);

  return 0.0;
}

void tabulate_kernel(void)
{
  int i;

  float lu,lumin,lumax,dlu,u;

  lumin=-6.;
  lumax=log10(2.1);
  dlu=(lumax-lumin)/(float)(NTAB-1);

  i=0;
  lu=lumin;
  
  while(i<NTAB)
    {
      u=pow(10.,lu);
      rh_table[i]=u;
      kernel_table[i]=cubic_spline_kernel(u);
      lu+=dlu;
      i++;
    }

}

float get_kernel_value(float u)
{
  float kval;
  
  int nmin=0,nmid,nmax=NTAB-1;

  if(u==rh_table[nmin])
    return kernel_table[nmin];

  if(u==rh_table[nmax])
    return kernel_table[nmax];  
  
  while((rh_table[nmin]-u)*(rh_table[nmax]-u)<0)
    {
      nmid=(nmin+nmax)/2;
      if((rh_table[nmid]-u)*(rh_table[nmax]-u)<0)
	nmin=nmid;
      else
	nmax=nmid;
      if(nmax-nmid<=1) break;
    }
  
  if(u==rh_table[nmin])
    return kernel_table[nmin];

  if(u==rh_table[nmax])
    return kernel_table[nmax];  

  kval=(rh_table[nmax]-u)/(rh_table[nmax]-rh_table[nmin])*kernel_table[nmin]+
    (u-rh_table[nmin])/(rh_table[nmax]-rh_table[nmin])*kernel_table[nmax];

  return kval;
  
}
    
float cumulative_cubic_spline_interpolant(float u)
{
  float sum=0.0;

  if(u<1.0)
    sum=u-0.5*u*u*u+(3./16.)*u*u*u*u;
  else if(u>=1. && u<2.)
    sum=0.6875+2*u-1.5*u*u+0.5*u*u*u-(1./16)*u*u*u*u-0.9375;
  else
    sum=0.75;

  sum/=M_PI;
      
  return sum;
}
  
void tabulate_integral(void)
{
  int i;

  float lu,lumin,lumax,dlu,u;

  lumin=-2.;
  lumax=log10(30.);
  dlu=(lumax-lumin)/(float)(NTAB-1);

  i=0;
  lu=lumin;
  
  while(i<NTAB)
    {
      u=pow(10.,lu);
      rh_table[i]=u;
      integral_table[i]=cumulative_cubic_spline_interpolant(u);
      lu+=dlu;
      i++;
    }

}

float get_integral_value(float u)
{
  float kval;
  
  int nmin=0,nmid,nmax=NTAB-1;

  if(u==rh_table[nmin])
    return integral_table[nmin];

  if(u==rh_table[nmax])
    return integral_table[nmax];  
  
  while((rh_table[nmin]-u)*(rh_table[nmax]-u)<0)
    {
      nmid=(nmin+nmax)/2;
      if((rh_table[nmid]-u)*(rh_table[nmax]-u)<0)
	nmin=nmid;
      else
	nmax=nmid;
      if(nmax-nmid<=1) break;
    }
  
  if(u==rh_table[nmin])
    return integral_table[nmin];

  if(u==rh_table[nmax])
    return integral_table[nmax];  

  kval=(rh_table[nmax]-u)/(rh_table[nmax]-rh_table[nmin])*integral_table[nmin]+
    (u-rh_table[nmin])/(rh_table[nmax]-rh_table[nmin])*integral_table[nmax];

  return kval;
  
}


