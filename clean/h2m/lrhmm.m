function [A_, mu_, Sigma_, logl] = lrhmm (X, st, A, mu, Sigma, n_it)

%lrhmm     Performs multiple iterations of the EM algorithm for a left-right
%	   model (with multiple training sequences).
%	   This is just a call to the lower-level functions.
%          Use: [A_,mu_,Sigma_,logl] = lrhmm(X,st,A,mu,Sigma,n_it)

% H2M Toolbox, Version 2.0
% Olivier Cappé, 07/02/97 - 26/02/98
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Turn verbose mode off
QUIET = 1;
% Check input data and dimensions
error(nargchk(6, 6, nargin));
% Check that the initialization value are correct
[N, p, DIAG_COV] = hmm_chk(A, [1 zeros(1, length(A(1,:))-1)], mu, Sigma);
% Use supplied values as initialization
A_ = A;
mu_ = mu;
Sigma_ = Sigma;

logl = zeros(1,n_it);
for i = 1:n_it
  [A_, logl(i), gamma] = hmm_mest(X, st, A_, mu_, Sigma_, QUIET);
  [mu_,Sigma_] = mix_par(X, gamma, DIAG_COV, QUIET);
  fprintf(1, 'Iteration %d:\t%.3f\n', i, logl(i));
end
