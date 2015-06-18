function [A_, pi0_, mu_, Sigma_, logl] = hmm (X, A, pi0, mu, Sigma, n_it)

%hmm       Performs multiple iterations of the EM algorithm.
%	   This is just a call to the lower-level functions.
%          Use: [A_,pi0_,mu_,Sigma_,logl] = hmm(X,A,pi0,mu,Sigma,n_it)

% H2M Toolbox, Version 2.0
% Olivier Cappé, 07/02/97 - 14/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Turn verbose mode off
QUIET = 1;
% Check input data and dimensions
error(nargchk(6, 6, nargin));
% Check that the initialization value are correct
[N, p, DIAG_COV] = hmm_chk(A, pi0, mu, Sigma);
% Check size of training data
[T,v] = size(X);
if (v ~= p)
  error('Size of training dat is incorrect.');
end 
% Use supplied values as initialization
A_ = A;
pi0_ = pi0;
mu_ = mu;
Sigma_ = Sigma;

logl = zeros(1,n_it);
for i = 1:n_it
  [alpha, beta, logl_tmp, dens] = hmm_fb(X, A_, pi0_, mu_, Sigma_, QUIET);
  logl(i) = logl_tmp;
  [A_, pi0_] = hmm_tran(alpha, beta, dens, A_, pi0_, QUIET);
  [mu_, Sigma_] = hmm_dens(X, alpha, beta, DIAG_COV, QUIET);
  fprintf(1, 'Iteration %d:\t%.3f\n', i, logl(i));
end
