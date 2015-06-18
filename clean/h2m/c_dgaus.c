/************************************************************************/
/* c_dgaus      Mex-file for computing a set of multivariate normal     */
/*              density values in the case of diagonal covariance       */
/*              matrices.                                               */
/*                                                                      */
/*      MATLAB 5 VERSION (OBSOLETE MATLAB 4 CALLS ARE INDICATED AS      */
/*      COMMENTS WHERE APPROPRIATE - 2 lines)                           */
/*                                                                      */
/* H2M Toolbox, Version 2.0                                             */
/* Olivier Cappé, 12/03/96 - 09/06/99                                   */
/* ENST Dpt. Signal / CNRS URA 820, Paris                               */
/************************************************************************/
# include <math.h>
# include "mex.h"

void compute_diag_gaus (double *X, double *mu, double *Sigma, double *dens,
  unsigned int N, unsigned int T, unsigned int p);

/************************************************************************/
/* MATLAB standard interface function.                                  */
/************************************************************************/
/* REPLACE BY THE FOLLOWING LINE IF YOU ARE STIL USING MATLAB V 4          */
/* void mexFunction (int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[]) { */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *X;
  double *mu;
  double *Sigma;
  double *dens;
  unsigned int N, T, p, test;

  if (nrhs != 3) {
    mexErrMsgTxt("Three input argument are required.");
  }
  /* Get dimensions */
  T = mxGetM(prhs[0]);
  p = mxGetN(prhs[0]);
  N = mxGetM(prhs[1]);
  test = mxGetN(prhs[1]);
  if (test != p) {
    mexErrMsgTxt("Wrong dimension for mean vectors.");
  }
  test = mxGetM(prhs[2]);
  if (test != N) {
    mexErrMsgTxt("Wrong number of diagonal covariances.");
  }
  test = mxGetN(prhs[1]);
  if (test != p) {
    mexErrMsgTxt("Wrong dimension for diagonal covariances.");
  }
  /* Pointers to array */
  X = mxGetPr(prhs[0]);
  mu = mxGetPr(prhs[1]);
  Sigma = mxGetPr(prhs[2]);
  /* Outupt values */
  /* REPLACE BY THE FOLLOWING LINE IF YOU ARE STIL USING MATLAB V 4        */
  /*  plhs[0] = mxCreateFull(T, N, REAL);                                  */
  plhs[0] = mxCreateDoubleMatrix(T, N, mxREAL);
  dens = mxGetPr(plhs[0]);
  /* Density computation */
  compute_diag_gaus (X, mu, Sigma, dens, N, T, p);
}

/************************************************************************/
/* Computational routine.                                               */
/************************************************************************/
void compute_diag_gaus (double *X, double *mu, double *Sigma, double *dens,
  unsigned int N, unsigned int T, unsigned int p) {
  unsigned int n, ln, t, lt, i, li_X, li_mu;
  double norm;
  double det;

  /* Normalizing factor */
  norm = pow((double) 6.28318530717959, - (double) p * (double) 0.5);
  /* Loop on mixture index */
  for (n=0, ln=0; n<N; n++, ln+=T) {
    /* Compute sqrt of det */
    det = (double) 1.0;
    for (i=0, li_mu=n; i<p; i++, li_mu+=N) {
      det *= Sigma[li_mu];
    }
    det = sqrt(det);
    /* Loop on time index */
    for (t=0, lt=ln; t<T; t++, lt++) {
      dens[lt] = (double) 0.0;
      for (i=0, li_X=t, li_mu=n; i<p; i++, li_X+=T, li_mu+=N) {
        dens[lt] += (X[li_X] - mu[li_mu]) * (X[li_X] - mu[li_mu])
          / Sigma[li_mu];
      }
      dens[lt] = norm * exp(- (double) 0.5 * dens[lt]) / det;
    }
  }
}
