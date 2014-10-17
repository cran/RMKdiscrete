#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

//Array of values separated by amount that should be greater than double eps:
const double pcutz[21] = {2, 3e-16, 3e-32, 3e-48, 3e-64, 3e-80, 3e-96, 3e-112, 3e-128, 3e-144, 3e-160, 3e-176, 3e-192,
  3e-208, 3e-224, 3e-240, 3e-256, 3e-272, 3e-288, 3e-304, 3e-308};

//*phold_ra points to a 21-length array that holds partial sums of probabilities
//newp is a new probability term to be added to the sum
void carefulprobsum(double newp, double *phold_ra, int add_carefully){
  int i, stopflag=0;
  //if adding carefully, we find the last element of pcutz that's greater than newp and add newp to the 
  //corresponding element of the phold array
  if(add_carefully==1){
    for(i=20;i>=0 && stopflag==0;i--){
      if(newp<pcutz[i]){
        phold_ra[i] += newp;
        stopflag = 1;
    }}
    //clean up the phold array so that each element is smaller than the corresponding value of pcutz
    for(i=20;i>0;i--){
      if(phold_ra[i]>pcutz[i]){
        phold_ra[i-1] += phold_ra[i];
        phold_ra[i] = 0;
    }}
  }
  //if not adding carefully, use only the first two elements of pcutz and the phold array:
  else{
    for(i=1;i>=0 && stopflag==0;i--){
      if(newp<pcutz[i]){
        phold_ra[i] += newp;
        stopflag = 1;
    }}
    if(phold_ra[1]>pcutz[1]){
      phold_ra[0] += phold_ra[1];
      phold_ra[1] = 0;
  }}
}

//Function to sum across the elements of phold array, from least to greatest, to return a scalar sum of 
//probabilities:
double carefulprobsum_fin(double *phold_ra, int add_carefully)
{
  int i;
  double out=0;
  //if adding carefully, need to use all 21 elements; otherwise use only first two:
  if(add_carefully==1){
    for(i=20;i>=0;i--){out += phold_ra[i];}
    return(out);
  }
  else{
    for(i=1;i>=0;i--){out += phold_ra[i];}
    return(out);
  }
}
