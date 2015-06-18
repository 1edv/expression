function [A_, logl, gamma] = hmm_mest (X, st, A, mu, Sigma, QUIET)

%hmm_mest  Reestimates transition parameters for multiple observation sequences
%          (left-to-right HMM).
%	Use: [A_,logl,gamma] = hmm_mest(X,st,A,mu,Sigma)
%	gamma can be used to reestimate the gaussian parameters.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 22/03/96 - 04/03/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Check input data and dimensions
error(nargchk(5, 6, nargin));
if (nargin < 6)
  QUIET = 0;
end
N = length(A(1,:));
pi0 = [1 zeros(1,N-1)];
[N,p,DIAG_COV] = hmm_chk(A, pi0, mu, Sigma);
[KT, nc] = size(X);
if (nc ~= p)
  error('Dimension of the observation vectors is incorrect.');
end
K = length(st);
if ((st(1) ~= 1) | (st(K) >= KT))
  error('Check the limits of the observation sequences.');
end

if (~QUIET)
  fprintf(1, 'Reestimating transition parameters..'); time = cputime;
end
% Force A to be sparse
A = sparse(A);
A_ = sparse(zeros(N));
gamma = zeros(KT, N);
den = 0;
logl = 0;
for i=1:K
  if (~QUIET)
    if (rem(i,5) == 1)
      fprintf(1, '.');
    end
  end
  % Limits of the observation sequence
  first = st(i);
  if (i == K)
    last = KT;
  else
    last = st(i+1)-1;
  end
  T = last - first + 1;
  [alpha, beta, logl_tmp, dens] = hmm_fb (X(first:last,:), A, pi0, mu, Sigma, 1);
  % Compute log-likelihood
  logl = logl + logl_tmp;
  % Reestimate transition parameters
  % Numerator
  for j=1:T-1
    A_i = alpha(j,:)' * (dens(j+1,:) .* beta(j+1,:)) .* A;
    A_i = A_i /sum(sum(A_i));
    A_ = A_ + A_i;
  end
  % Denominator
  % Memorize state probabilities
  gamma(first:last,:) = alpha .* beta;
  gamma(first:last,:) = gamma(first:last,:) ./ (sum(gamma(first:last,:)')' * ones(1,N));
  den = den + sum(gamma(first:first+T-2,:));
end
% Normalize transition matrix
A_ = A_ ./ (den' * ones(1,N));
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n', time);
end
