function [alpha, beta, logl, dens] = hmm_fb (X, A, pi0, mu, Sigma, QUIET)

%hmm_fb    Implements the forward-backward recursions (with scaling).
%	Use: [alpha,beta,logl,dens] = hmm_fb(X,A,pi0,mu,Sigma) where
%	alpha and beta are the forward and backward variables and logl
%	is the log-likelihood). dens contains the values of all
%	gaussian densities for each time index (usefull for the
%	estimation of the transition probabilities).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 11/03/96 - 14/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Input args.
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
dens = gauseval(X, mu, Sigma, QUIET);

% Forward pass
if (~QUIET)
  fprintf(1, 'Forward pass...'); time = cputime;
end
scale = [1 ; zeros(T-1, 1)];
alpha = zeros(T, N);
alpha(1,:) = pi0.*dens(1,:);
for i=2:T
  alpha(i,:) = (alpha(i-1,:) * A) .* dens(i,:);
  % Systematic scaling
  scale(i) = sum(alpha(i,:));
  alpha(i,:) = alpha(i,:) / scale(i);
end
% This is not computationnaly efficient but log(prod()) won't work here !
logl = sum(log(scale));
if (~QUIET)
  time = cputime - time;
  fprintf(1, ' (%.2f s)\n', time);
end

% Backward pass
if (~QUIET)
  fprintf(1, 'Backward pass...'); time = cputime;
end
% Scale the backward variable with the forward scale factors (this ensures
% that the reestimation of the transition matrix is correct)
beta = zeros(T, N);
beta(T, :) = ones(1, N);
for i=(T-1):-1:1
  beta(i, :) = (beta(i+1,:).* dens(i+1,:)) * A';
  % Apply scaling
  beta(i,:) = beta(i,:) / scale(i);
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n',time);
end
