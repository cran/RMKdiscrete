//Declare functions that are used in more than one file of C source:

void carefulprobsum(double newp, double *phold_ra, int add_carefully);

double carefulprobsum_fin(double *phold_ra, int add_carefully);

double do_LGP_findmax(double theta, double lambda);

double do_dLGP(double x, double theta, double lambda, double nc, int give_log);

double do_dLGP_withmax(double x, double theta, double lambda, double nc, int give_log, double max);

double do_LGP_getnc(double nctol, double theta, double lambda, int add_carefully);
