#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "00000000.h"


//Function to compute bivariate LGP PMF for one pair (x,y):
double do_dbiLGP(double x, double y, double theta0, double theta1, double theta2, double lambda0, double lambda1, 
  double lambda2, double nc0, double nc1, double nc2, int give_log, int add_carefully)
{
  double out, phold;
  if(theta0==0){ //Stochastic independence
    out = do_dLGP(x,theta1,lambda1,nc1,1) + do_dLGP(y,theta2,lambda2,nc2,1);
    return( give_log==1 ? out : exp(out) );
  }
  out = phold = 0;
  double u;
  double umax = fmin2(fmin2(x,y), do_LGP_findmax(theta0,lambda0));
  double max_xp = do_LGP_findmax(theta1,lambda1);
  double max_yp = do_LGP_findmax(theta2,lambda2);
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(u=0;u<=umax;u++){
    /*From a numerical stability standpoint, is doing this multiplication by adding logs and exponentiating the sum
    really the smartest way to do it?*/
    phold = exp(do_dLGP_withmax(x-u,theta1,lambda1,nc1,1,max_xp) + do_dLGP_withmax(y-u,theta2,lambda2,nc2,1,max_yp) + 
      do_dLGP_withmax(u,theta0,lambda0,nc0,1,umax));
    carefulprobsum(phold, parray, add_carefully);
    R_CheckUserInterrupt();
  }
  out = carefulprobsum_fin(parray,add_carefully);
  out = ((give_log==1) ? log(out) : out);
  return(out);
}

//Function that manages inputs from frontend, and invokes do_dbiLGP() while looping through:
void call_dbiLGP(double *x, double *y, double *theta0, double *theta1, double *theta2, double *lambda0, double *lambda1,
  double *lambda2, double *nc0, double *nc1, double *nc2, int *give_log, int *add_carefully, int *Cnout, double *Cout)
{
    int i;
    for(i=0;i<*Cnout;i++){ 
      Cout[i] = do_dbiLGP(x[i],y[i],theta0[i],theta1[i],theta2[i],lambda0[i],lambda1[i],
        lambda2[i],nc0[i],nc1[i],nc2[i],*give_log,*add_carefully);
      R_CheckUserInterrupt();
    }
}

//Function to calculate PMF of sum of two independent LGP r.v.s:
double do_dLGP_convolution(double x, double theta0, double theta1, double lambda0, double lambda1, double nc0, 
double nc1, int add_carefully) 
{
  double out = 0;
  if(lambda0==lambda1 && lambda0>=0){ //Addition rule applies if lambdas are both equal to some non-neg value
    //double nc = do_LGP_getnc(1e-14,theta0+theta1,lambda0,add_carefully);
    return(do_dLGP(x,theta0+theta1,lambda0,1.0,0L));
  }
  double u = 0, phold = 0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double umax = fmin2(x, do_LGP_findmax(theta0,lambda0));
  for(u=0;u<=umax;u++){
    phold = exp(do_dLGP(x-u,theta1,lambda1,nc1,1)+do_dLGP_withmax(u,theta0,lambda0,nc0,1,umax));
    carefulprobsum(phold, parray, add_carefully);
  }
  out = carefulprobsum_fin(parray,add_carefully);
  return(out);
}

/*Function to numerically calculate means, variances, and covariance of log-transformed bivariate LGP distribution
(looping over inputs is done in R):*/
void call_biLGP_logMV(double *theta0, double *theta1, double *theta2, 
  double *lambda0, double *lambda1, double *lambda2,
  double *nc0, double *nc1, double *nc2,
  double *const_add, double *tol, int *add_carefully,
  double *EX, double *EY, double *EX2, double *EY2, double *EXY)
{ //What are called Y1 and Y2 in the documentation are here identified as X and Y...
    double nexterm=0, oldterm=0;
    int xmodeflag=0;
    int xstopflag=0;
    double i=0, j=0, x, y;
    if(*lambda0==1 || *lambda1==1){//Moments are not finite when lambda=1.
      *EX = R_PosInf;
      *EX2 = R_PosInf;
    }
    else{
      for(i=0;xstopflag==0;i++){
        nexterm = do_dLGP_convolution(i,*theta0,*theta1,*lambda0,*lambda1,*nc0,*nc1,*add_carefully);
        if(nexterm < oldterm) xmodeflag = 1; //if mode of x's distribution has been crossed
        *EX += nexterm * log(i + *const_add);
        *EX2 += nexterm * R_pow_di(log(i + *const_add),2L);
        //stop if second raw moment has converged:
        if(nexterm * R_pow_di(log(i + *const_add),2L) < *tol && xmodeflag==1) xstopflag=1; 
        //if(nexterm==0) xstopflag=1;
        oldterm = nexterm;
      }
    }
    R_CheckUserInterrupt();
    //Now do for y as was done for x, unless they have the same marginal distributions:
    if( *theta1==*theta2 && *lambda1==*lambda2 ){
      *EY = *EX;
      *EY2 = *EX2;
      j = i;
    }
    else{
      int ymodeflag=0, ystopflag=0;
      oldterm=0;
      if(*lambda0==1 || *lambda2==1){
        *EY = R_PosInf;
        *EY2 = R_PosInf;
      }
      else{
        for(j=0;ystopflag==0;j++){
          nexterm = do_dLGP_convolution(j,*theta0,*theta2,*lambda0,*lambda2,*nc0,*nc2,*add_carefully);
          if(nexterm < oldterm) ymodeflag = 1;
          *EY += nexterm * log(j + *const_add);
          *EY2 += nexterm * R_pow_di(log(j + *const_add),2L);
          if(nexterm * R_pow_di(log(j + *const_add),2L) < *tol && ymodeflag==1) ystopflag=1;
          //if(nexterm==0) ystopflag=1;
          oldterm = nexterm;
    }}}
    R_CheckUserInterrupt();
    if( (*lambda0==1 || *lambda1==1) || *lambda2==1 ){*EXY = R_PosInf;}  //Again, check for infinite moments
    //compute E(XY), going as high on X and Y as was done for their marginal second moments:
    else{
      for(x=0;x<=i;x++){
        for(y=0;y<=j;y++){
          *EXY += do_dbiLGP(x,y,*theta0,*theta1,*theta2,*lambda0,*lambda1,*lambda2,*nc0,*nc1,*nc2,0,*add_carefully) * 
            log(x + *const_add) * log(y + *const_add);
        }
      R_CheckUserInterrupt();
      }
}}
