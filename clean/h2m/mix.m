function [w_, mu_, Sigma_, logl] = mix (X, w, mu, Sigma, n_it)

%mix       Performs multiple iterations of the EM algorithm for a mixture model.
%	   This is just a call to the lower-level functions.
%          Use: [w_,mu_,Sigma_,logl] = mix (X,w,mu,Sigma,n_it)

% H2M Toolbox, Version 2.0
% Olivier Cappé, 07/02/97 - 04/03/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Turn verbose mode off
QUIET = 1;
% Check input data and dimensions
error(nargchk(5, 5, nargin));
% Check that the initialization value are correct
[N, p, DIAG_COV] = mix_chk(w, mu, Sigma);
% Use supplied values as initialization
w_ = w;
mu_ = mu;
Sigma_ = Sigma;

logl = zeros(1,n_it);
for i = 1:n_it
  [gamma, logl(i)] = mix_post(X, w_, mu_, Sigma_, QUIET);
  [mu_, Sigma_, w_] = mix_par(X, gamma, DIAG_COV, QUIET);
  fprintf(1, 'Iteration %d:\t%.3f\n', i, logl(i));
end
