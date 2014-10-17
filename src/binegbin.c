#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "00000000.h"


//Function to compute bivariate negbin PMF for one pair (x,y):
double do_dbinegbin(double x, double y, double nu0, double nu1, double nu2, double p0, double p1, 
  double p2, int give_log, int add_carefully){
  double out, phold;
  if(nu0==0){
    out = dnbinom(x,nu1,p1,1) + dnbinom(y,nu2,p2,1);
    return( give_log==1 ? out : exp(out) );
  }
  out = phold = 0;
  double u;
  double umax = fmin2(x,y);
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(u=0;u<=umax;u++){
    phold = exp(dnbinom(x-u,nu1,p1,1)+dnbinom(y-u,nu2,p2,1)+dnbinom(u,nu0,p0,1));
    carefulprobsum(phold, parray, add_carefully);
    R_CheckUserInterrupt();
  }
  out = carefulprobsum_fin(parray,add_carefully);
  out = ((give_log==1) ? log(out) : out);
  return(out);
}

//Function that manages inputs from frontend, and invokes do_dbinegbin() while looping through:
void call_dbinegbin(double *x, double *y, double *nu0, double *nu1, double *nu2, double *p0, double *p1, 
  double *p2, int *give_log, int *add_carefully, int *Cnout, double *Cout){
    int i;
    for(i=0;i<*Cnout;i++){ 
      Cout[i] = do_dbinegbin(x[i],y[i],nu0[i],nu1[i],nu2[i],p0[i],p1[i],p2[i],*give_log,*add_carefully);
      R_CheckUserInterrupt();
    }
}

double do_dnegbin_convolution(double x, double nu0, double nu1, double p0, double p1, int add_carefully){
  double out = 0;
  if(p0==p1){
    return(dnbinom(x,nu0+nu1,p0,0));
  }
  double u = 0, phold = 0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(u=0;u<=x;u++){
    phold = exp(dnbinom(x-u,nu1,p1,1L)+dnbinom(u,nu0,p0,1L));
    carefulprobsum(phold, parray, add_carefully);
  }
  out = carefulprobsum_fin(parray,add_carefully);
  return(out);
}

void call_binegbin_logMV(double *nu0, double *nu1, double *nu2, 
  double *p0, double *p1, double *p2,
  double *const_add, double *tol, int *add_carefully,
  double *EX, double *EY, double *EX2, double *EY2, double *EXY){
    double nexterm=0, oldterm=0;
    int xmodeflag=0;
    int xstopflag=0;
    double i=0, j=0, x, y;
    for(i=0;xstopflag==0;i++){
      nexterm = do_dnegbin_convolution(i,*nu0,*nu1,*p0,*p1,*add_carefully);
      if(nexterm < oldterm) xmodeflag = 1;
      *EX += nexterm * log(i + *const_add);
      *EX2 += nexterm * R_pow_di(log(i + *const_add),2);
      if(nexterm * R_pow_di(log(i + *const_add),2) < *tol && xmodeflag==1) xstopflag=1;
      //if(nexterm==0) xstopflag=1;
      oldterm = nexterm;
    }
    R_CheckUserInterrupt();
    //Now do for y as was done for x, unless they have the same marginal distributions:
    if( *nu1==*nu2 && *p1==*p2 ){
      *EY = *EX;
      *EY2 = *EX2;
      j = i;
    }
    else{
      int ymodeflag=0, ystopflag=0;
      oldterm=0;
      for(j=0;ystopflag==0;j++){
        nexterm = do_dnegbin_convolution(j,*nu0,*nu2,*p0,*p2,*add_carefully);
        if(nexterm < oldterm) ymodeflag = 1;
        *EY += nexterm * log(j + *const_add);
        *EY2 += nexterm * R_pow_di(log(j + *const_add),2);
        if(nexterm * R_pow_di(log(j + *const_add),2) < *tol && ymodeflag==1) ystopflag=1;
        //if(nexterm==0) ystopflag=1;
        oldterm = nexterm;
      }}
    R_CheckUserInterrupt();
    for(x=0;x<=i;x++){
      for(y=0;y<=j;y++){
        *EXY += do_dbinegbin(x,y,*nu0,*nu1,*nu2,*p0,*p1,*p2,0,*add_carefully) * 
          log(x + *const_add) * log(y + *const_add);
        }
      R_CheckUserInterrupt();
    }
}
