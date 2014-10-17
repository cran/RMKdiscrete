#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "00000000.h"

/*Functions designated with "do_" are workhorses that get called by other C functions. 
Functions designated with "call_" are what the frontend
interfaces with; they loop through inputs and do their thing. */

double do_LGP_findmax(double theta, double lambda)
{
  if(theta<0 || fabs(lambda)>1){return(R_NaN);} //should always be checked in R, but...
  if(theta==0){return(0);}
  if(lambda>=0){return(R_PosInf);}
  double max = 1;
  if(-theta/lambda <= 1){max = 0;}
  if(-theta/lambda > 1){
      max = ftrunc(-theta/lambda);
      if(ftrunc(-theta/lambda) == -theta/lambda) max-- ; //Has to be NEXT LOWEST integer
  }
  return(max);
}

void call_LGP_findmax(double *theta, double *lambda, int *Cnout, double *Cout)
{
  int i;
  for(i=0;i<*Cnout;i++){ //i<*Cnout
    Cout[i] = do_LGP_findmax(theta[i],lambda[i]); 
    R_CheckUserInterrupt();
  } 
}

double do_dLGP(double x, double theta, double lambda, double nc, int give_log)
{
  double out;
  double max;
  if(theta==0 && x==0) return( (give_log==1) ? 0 : 1 );
  if(theta==0 && x!=0) return( (give_log==1) ? R_NegInf : 0 );
  if(lambda==0) return( dpois(x,theta,give_log) );
  if(lambda < 0){
    max = do_LGP_findmax(theta,lambda);
    if(x>max) return( (give_log==1) ? R_NegInf : 0 );
  }
  /*It's possible that the following method of computing the PMF is not the smartest--see Catherine Loader's 
  (2000) paper on the binomial and Poisson*/
  out = log(theta) + (x-1)*log(theta+(lambda * x)) - theta - (lambda * x) - lgammafn(x+1) - log(nc);
  return( (give_log==1) ? out : exp(out) );
}

//Function to evaluate PMF in C when upper support limit has already been found for a different reason:
//(TOmaybeDO: combine this with do_dLGP())
double do_dLGP_withmax(double x, double theta, double lambda, double nc, int give_log, double max)
{
  double out;
  if(theta==0 && x==0) return( (give_log==1) ? 0 : 1 );
  if(theta==0 && x!=0) return( (give_log==1) ? R_NegInf : 0 );
  if(lambda==0) return( dpois(x,theta,give_log) );
  if(x>max) return( (give_log==1) ? R_NegInf : 0 );
  out = log(theta) + (x-1)*log(theta+(lambda * x)) - theta - (lambda * x) - lgammafn(x+1) - log(nc);
  return( (give_log==1) ? out : exp(out) );
}

void call_dLGP(double *x, double *theta, double *lambda, double *nc, int *give_log, int *Cnout, double *Cout)
{
  int i;
  for(i=0;i<*Cnout;i++){
    Cout[i] = do_dLGP(x[i],theta[i],lambda[i],nc[i],*give_log); 
    R_CheckUserInterrupt();
  }
}

double do_LGP_getnc(double nctol, double theta, double lambda, int add_carefully)
{
  if(lambda >= 0) return(1); //no need to numerically normalize the probabilities for non-negative lambda
  double max = do_LGP_findmax(theta,lambda);
  if(max==0) return(do_dLGP_withmax(0,theta,lambda,1,0,0));
  nctol = (max>200000) ? nctol : 0; /*nctol is ignored unless max is big, when it might be needed to make the 
                                    loop stop within a reasonable amount of time*/
  double z, nc = 0, holder=0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(z=0;z<=max && fabs(1-nc)>nctol;z++){ //Add probabilities until max is reached or nc is suff. close to 1
    holder = do_dLGP_withmax(z,theta,lambda,1,0,max);
    carefulprobsum(holder, parray, add_carefully);
    R_CheckUserInterrupt();
  }
  nc = carefulprobsum_fin(parray, add_carefully);
  return(nc);
}

void call_LGP_getnc(double *nctol, double *theta, double *lambda, int *Cnout, double *Cout, int *add_carefully)
{
  int i;
  for(i=0;i<*Cnout;i++){
    Cout[i] = do_LGP_getnc(*nctol,theta[i],lambda[i],*add_carefully);
    R_CheckUserInterrupt();
  }
}

//Function for LGP summary when lambda is negative:
void call_sLGP_neglam(double *theta, double *lambda, double *nc, int *Cnout, double *mu1, double *med, double *mod,
double *mu2, double *mu3, double *mu4, int *add_carefully)
{
  //Initialize:
  double nexterm=0, lognexterm=0, logoldterm=0, max=0, x=0, logx=0, EX=0, EX2=0, EX3=0, EX4=0, EX4_term=0;
  int i=0, k=0, modeflag=0, medianflag=0, stopflag=0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(i=0;i<*Cnout;i++){ //Loop across inputs passed from frontend
    max = do_LGP_findmax(theta[i],lambda[i]);
    modeflag = 0, medianflag = 0, stopflag=0;
    //Rules for mode-finding, from Consul (1989):
    if( theta[i] * exp(-1*lambda[i]) < 1 ){
      mod[i] = 0;
      modeflag = 1; //modeflag=1 means that the mode has been found
    }
    if( theta[i] * exp(-1*lambda[i]) == 1 ){
      mod[i] = 0.5;
      modeflag = 1;
    }
    logoldterm=R_NegInf, EX=0, EX2=0, EX3=0, EX4=0, EX4_term=0;
    for(k=0;k<21;k++){parray[k] = 0;} //Need to re-zero parray every iteration.
    for(x=0;x<=max && stopflag==0;x++){
      logx = log(x);
      lognexterm = do_dLGP_withmax(x,theta[i],lambda[i],nc[i],1,max);
      nexterm = exp(lognexterm);
      //If mode not yet found and probabilities are now decreasing, the previous x must have been the mode:
      if( lognexterm < logoldterm && modeflag==0){
        mod[i] = x-1;
        modeflag = 1;
      }
      //if median not yet found, cumulative probability is computed and stored
      if(medianflag==0){
        carefulprobsum(nexterm, parray, *add_carefully);
        //median found once cumu. prob. >= 0.5
        if( carefulprobsum_fin(parray, *add_carefully) >= 0.5 ){
          med[i] = x;
          medianflag = 1;
        }
      }
      //TODO: Improve numerical precision of the 3rd & 4th moments.
      EX += exp(lognexterm + logx);
      EX2 += exp(lognexterm + 2*logx);
      EX3 += exp(lognexterm + 3*logx);
      EX4_term = exp(lognexterm + 4*logx);
      EX4 += EX4_term;
      logoldterm = lognexterm;
      R_CheckUserInterrupt();
      //if median and mode found and new terms added to 4th moment are underflowing to zero, then stop:
      if(modeflag==1 && medianflag==1 && EX4_term==0){stopflag = 1;}
    }
    //if mode still not found, presumably it is the upper support limit:
    if(modeflag==0){mod[i] = max;}
    //Compute central moments from raw moments:
    mu1[i] = EX;
    mu2[i] = EX2 - (EX * EX);
    mu3[i] = EX3 - (3 * EX2 * EX) + (2 * EX * EX * EX);
    //printf("%f \n %f \n %f \n %f \n",EX,EX2,EX3,EX4);
    mu4[i] = 
      EX4 - (EX3 * EX) +
      (EX2 * EX * EX) + (EX2 * EX * EX) + (EX2 * EX * EX) - 
      (EX3 * EX) + (EX2 * EX * EX) - (EX * EX * EX * EX) -
      (EX3 * EX) + (EX2 * EX *EX) - (EX * EX * EX * EX) - 
      (EX3 * EX) + (EX2 * EX *EX) - (EX * EX * EX * EX);
  }
}

//Find the mode, for use with lambda >= 0.
void call_LGP_findmode(double *theta, double *lambda, double *nc, double *start, int *Cnout, double *Cout,
  int *failflag)
{
  //Initialize:
  double nexterm=0, oldterm=0, x=0, max=0;
  int i=0, stopflag=0;
  for(i=0;i<*Cnout;i++){ //Loop across inputs from frontend
    max = do_LGP_findmax(theta[i],lambda[i]);
    if(max==0){
      Cout[i] = 0;
      continue;
    }
    if( theta[i] * exp(-1*lambda[i]) < 1 ){
      Cout[i] = 0;
      continue;
    }
    if( theta[i] * exp(-1*lambda[i]) == 1 ){
      Cout[i] = 0.5;
      continue;
    }
    stopflag = 0;
    x = start[i]; //start is passed from frontend; is lower bound on mode, per Consul (1989)
    oldterm = do_dLGP_withmax(x,theta[i],lambda[i],nc[i],1,max);
    //Find probabilities of increasing x until they start to decrease, or max is reached:
    while(stopflag==0){
      x++;
      nexterm = do_dLGP_withmax(x,theta[i],lambda[i],nc[i],1,max);
      if(nexterm < oldterm){
        Cout[i] = x-1;
        stopflag = 1;
      }
      if(x==max){
        if(stopflag==0){Cout[i] = x;}
        stopflag = 1;
      }
      oldterm = nexterm;
      R_CheckUserInterrupt();
    }
  }
}

//pLGP backend, general purpose:
void call_pLGP(double *q, double *theta, double *lambda, double *nc, int *lower_tail, int *Cnout, 
  double *Cout, int *failflag, double *i_fail, int *add_carefully)
{
    //Initialize:
    int i=0, modeflag=0, k=0;
    double oldterm=0, nexterm=0, j=0, max=0;
    //If lower_tail, we are adding probabilities to zero; otherwise, we end up subtracting them from 1:
    double plusorminus = ((*lower_tail==1) ? 1 : -1);
    double p =  ((*lower_tail==1) ? 0 : 1);
    double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for(i=0;i<*Cnout;i++){ //loop over inputs from frontend
      max = do_LGP_findmax(theta[i],lambda[i]);
      oldterm = 0;
      for(k=0;k<21;k++){parray[k] = 0;} //need to re-zero every iteration
      modeflag = 0;
      for(j=0;j<=q[i] && failflag[i]==0;j++){
        nexterm = do_dLGP_withmax(j,theta[i],lambda[i],nc[i],0,max);
        carefulprobsum(nexterm, parray, *add_carefully);
        /*Need to know if mode has been reached, because if it has and probabilities start underflowing to
        zero, we know they won't start increasing:*/
        if(nexterm < oldterm) modeflag = 1;
        if(nexterm==0 && modeflag==1){
          failflag[i] = 1;
          p = p + (plusorminus * carefulprobsum_fin(parray, *add_carefully));
          i_fail[i] = j; //So that user can be warned where PMF was computationally zero.
        }
        if(j==q[i]){ p = p + (plusorminus * carefulprobsum_fin(parray, *add_carefully)); }
        oldterm = nexterm;
        R_CheckUserInterrupt();
      }
      //Probabilities could conceivably fall slightly outside [0,1] due to numerical imprecision...
      if(p<0) p = 0;
      if(p>1) p = 1;
      Cout[i] = p;
}}


/*cumulative probabilities, for use when multiple quantiles are provided but only one value each for theta and 
lambda.  This function avoids redundant addition of probabilities.*/
void call_pLGP_lowertailsearch(double *q, double *theta, double *lambda, double *nc,
  int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
  {
    //^^^Note than q is a vector of quantiles, ordered from least to greatest.
    double start = -1, p = 0, nexterm=0, oldterm=0, j;
    *max = ((*max<0) ? R_PosInf : *max); //Avoids providing infinite values to .C() in R
    //stopflagi=1 means stop iterating across input quantiles; stopflagj=1 means stop for the current quantile:
    int modeflag=0, stopflagi=0, stopflagj=0, i=0;
    double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //loop across quantiles:
    for(i=0;i<*Cnout;i++){ 
      //if q[i] is the same as the previous q, or if no more iterating over q, use the previous q's result:
      if(stopflagi==1 || (i>0 && q[i]==q[i-1])) Cout[i] = Cout[i-1];
      else{
        stopflagj=0;
        j=start;
        //For q values less than zero or greater than max, the cumulative probability is easy:
        if(q[i] < 0){
          Cout[i] = 0;
          stopflagj = 1;
          continue;
        }
        if(q[i] >= *max){
          Cout[i] = 1;
          stopflagj = 1;
          stopflagi = 1;
          continue;
        }
        //Work upward from start until q[i] is reached or PMF underflows to zero past the mode:
        for(j=start+1;stopflagj==0 && j<=q[i];j++){
          nexterm = do_dLGP_withmax(j,*theta,*lambda,*nc,0,*max);
          carefulprobsum(nexterm, parray, *add_carefully);
          if(nexterm < oldterm) modeflag = 1;
          if(nexterm==0 && modeflag==1){ //if mode has been reached and PMF is returning zero...
            *failflag = 1;
            stopflagi = 1;
            stopflagj = 1;
            p = carefulprobsum_fin(parray, *add_carefully);
            *i_fail = j;
          }
          if(j==q[i]){
            p = carefulprobsum_fin(parray, *add_carefully);
          }
          oldterm = nexterm;
          R_CheckUserInterrupt();
        }
        if(p<0) p = 0;
        if(p>1){
          p = 1;
          stopflagi = 1;
        }
        Cout[i] = p;
        start = q[i];
  }
}}

/*Upper-tail probabilities, for use when multiple quantiles are provided but only one value each for theta and 
lambda.  This function avoids redundant addition of probabilities.*/
void call_pLGP_uppertailsearch(double *q, double *theta, double *lambda, double *nc,
  int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
{
  double start = -1, p = 1, nexterm=0, oldterm=0, j;
  *max = ((*max<0) ? R_PosInf : *max); //Avoids providing infinite values to .C()
  int modeflag=0, stopflagi=0, stopflagj=0, i=0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(i=0;i<*Cnout;i++){
    if(stopflagi==1 || (i>0 && q[i]==q[i-1])) Cout[i] = Cout[i-1];
    else{
      stopflagj=0;
      j=start;
      if(q[i] < 0){
        Cout[i] = 1;
        stopflagj = 1;
        continue;
      }
      if(q[i] >= *max){
        Cout[i] = 0;
        stopflagj = 1;
        stopflagi = 1;
        continue;
      }
      for(j=start+1;stopflagj==0 && j<=q[i];j++){
        nexterm = do_dLGP_withmax(j,*theta,*lambda,*nc,0,*max);
        carefulprobsum(nexterm, parray, *add_carefully);
        if(nexterm < oldterm) modeflag = 1;
        if(nexterm==0 && modeflag==1){
          *failflag = 1;
          stopflagi = 1;
          stopflagj = 1;
          p = 1 - carefulprobsum_fin(parray, *add_carefully);
          *i_fail = j;
        }
        if(j==q[i]){
          p = 1 - carefulprobsum_fin(parray, *add_carefully);
        }
        oldterm = nexterm;
        R_CheckUserInterrupt();
      }
      if(p<0) p = 0;
      if(p>1){
        p = 1;
        stopflagi = 1;
      }
      Cout[i] = p;
      start = q[i];
  }}
}

/*Upper-tail probabilities, for use when multiple quantiles are provided but only one value each for theta and 
lambda, lambda is negative, and the given quantiles are all in the upper half of the support.*/
void call_pLGP_uppertailsearch_neglam(double *q, double *theta, double *lambda, double *nc,
  int *Cnout, double *Cout, int *failflag, double *i_fail, int *add_carefully, double *max)
{
  *max = ((*max<0) ? R_PosInf : *max); //Avoids providing infinite values to .C()
  double start = *max, p = 1, nexterm=0, j;
  int stopflagi=0, stopflagj=0, i=0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(i=0;i<*Cnout;i++){
    if(stopflagi==1 || (i>0 && q[i]==q[i-1])) Cout[i] = Cout[i-1];
    else{
      stopflagj=0;
      j=start;
      if(q[i] < 0){
        Cout[i] = 1;
        stopflagj = 1;
        stopflagi = 1;
        continue;
      }
      if(q[i] >= *max){
        Cout[i] = 0;
        stopflagj = 1;
        continue;
      }
      for(j=start;stopflagj==0 && j>q[i];j--){
        nexterm = do_dLGP_withmax(j,*theta,*lambda,*nc,0,*max);
        carefulprobsum(nexterm, parray, *add_carefully);
        if(nexterm==0 && *failflag==0) *failflag = 1;
        if(nexterm>0 && *failflag==1){
          *i_fail = j;
          *failflag = 2;
        }
        if(j==q[i]+1){p = 1 - carefulprobsum_fin(parray, *add_carefully);}
        R_CheckUserInterrupt();
      }
      if(p>1) p = 1;
      if(p<0){
        p = 0;
        stopflagi = 1;
      }
      Cout[i] = p;
      start = q[i];
}}}

//qLGP, general-purpose
void call_qLGP(double *p, double *theta, double *lambda, double *nc, int *Cnout, double *Cout, 
  int *failflag, double *i_fail, double *pcumu, int *add_carefully){
    int i=0, k=0, modeflag=0, stopflag=0;
    double j=0, nexterm=0, oldterm=0, max=0;
    double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for(i=0;i<*Cnout;i++){
      modeflag=0, stopflag=0, j=0, nexterm=0, oldterm=0;
      max = do_LGP_findmax(theta[i],lambda[i]);
      for(k=0;k<21;k++){parray[k] = 0;}
      if(p[i]==0){
        Cout[i] = 0;
        stopflag=1;
      }
      if(p[i]==1){
        Cout[i] = max;
        stopflag = 1;
      }
      for(j=0;stopflag==0;j++){
        Cout[i] = j;
        nexterm = do_dLGP_withmax(j,theta[i],lambda[i],nc[i],0,max);
        if(nexterm < oldterm) modeflag = 1;
        carefulprobsum(nexterm, parray, *add_carefully);
        if( carefulprobsum_fin(parray, *add_carefully) >= p[i] ){ stopflag=1; }
        if(nexterm==0 && modeflag==1){
          failflag[i] = 1;
          stopflag = 1;
          i_fail[i] = j;
        }
        oldterm = nexterm;
        R_CheckUserInterrupt();
      }
    pcumu[i] = carefulprobsum_fin(parray, *add_carefully);
    }
}

/*Function for quantiles, when multiple cumu probabilities provided, but only one value each for theta and 
lambda.  Avoids redundant addition of probabilities.*/
void call_qLGP_pvec(double *p, double *theta, double *lambda, double *nc,
  int *Cnout, double *Cout, int *failflag, double *i_fail, double *pcumu, int *add_carefully, double *max)
{
  *max = ((*max<0) ? R_PosInf : *max);
  double j = -1, nexterm=0, oldterm=0;
  int i=0, stopflagi=0, stopflagj=0, modeflag=0;
  double parray[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(i=0;i<*Cnout;i++){
    if(stopflagi==1 || (i>0 && p[i]==p[i-1])){
      Cout[i] = Cout[i-1];
      continue;
    }
    if(p[i]==1){
      Cout[i] = *max;
      stopflagi = 1;
      continue;
    }
    if(p[i]==0){
      Cout[i] = 0;
      continue;
    }
    if(*pcumu >= p[i]){
      Cout[i] = j;
      continue;
    }
    stopflagj=0;
    while(stopflagj==0){
      j++;
      nexterm = do_dLGP_withmax(j,*theta,*lambda,*nc,0,*max);
      if(nexterm < oldterm) modeflag = 1;
      carefulprobsum(nexterm, parray, *add_carefully);
      if( carefulprobsum_fin(parray, *add_carefully) >= p[i] ){ stopflagj=1; }
      if(nexterm==0 && modeflag==1){
        *failflag = 1;
        stopflagj = 1;
        stopflagi = 1;
        *i_fail = j;
      }
      oldterm = nexterm;
      R_CheckUserInterrupt();
    }
    Cout[i] = j;
    *pcumu = carefulprobsum_fin(parray, *add_carefully);
  }
}
