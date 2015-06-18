function [logl, path, plogl, logdens] = hmm_vit (X, A, pi0, mu, Sigma, QUIET)

%hmm_vit   Computes the most likely sequence of states (Viterbi DP algorithm).
%	Use: [logl,path,plogl,dens] = hmm_vit(X,A,pi0,mu,Sigma). The
%	computation is peformed using partial log-likelihoods.

%	The mex-file 'c_dgaus' will be used in the diagonal case if it is
%	found in MATLAB's search path.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 05/07/96 - 17/06/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

error(nargchk(5, 6, nargin));
if (nargin < 6)
  QUIET = 0;
end
[N, p, DIAG_COV] = hmm_chk(A, pi0, mu, Sigma);
[T, Nc] = size(X);
if (Nc ~= p)
  error('Observation vectors have an incorrect dimension.');
end

% Compute density values
logdens = gauslogv(X, mu, Sigma, QUIET);

% Viterbi recursion
if (~QUIET)
  fprintf(1, 'Viterbi recursion...'); time = cputime;
end
% HMM probabilities in log
A = log(A + realmin);
pi0 = log(pi0 + realmin);
% Partial loglikelihood array and bactracking array
plogl = zeros(T, N);
bcktr = zeros(T-1, N);
plogl(1,:) = pi0 + logdens(1,:);
for i=2:T
  [plogl(i,:), bcktr(i-1,:)] = max((plogl(i-1,:)' * ones(1,N)) + A);
   plogl(i,:) = plogl(i,:) + logdens(i,:);
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n',time);
end

% Backtracking
if (~QUIET)
  fprintf(1, 'Backtracking...'); time = cputime;
end
path = zeros(T,1);
[logl, path(T)] = max(plogl(T,:));
for i=(T-1):-1:1
  path(i) = bcktr(i, path(i+1));
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n',time);
end
