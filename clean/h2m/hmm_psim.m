function s = hmm_psim (X, A, pi0, mu, Sigma, QUIET)

%hmm_psim  Generates a random sequence of conditional HMM states.
%	Use: s = hmm_psim(X,A,pi0,mu,Sigma) where s is a sequence of states
%	of the HMM drawn conditionnaly to the observed vectors X.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 11/07/97 - 31/07/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (~((exist('randindx') == 2) | (exist('randindx') == 3)))
  error('Function randindx is missing.');
end

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
alpha = zeros(T, N);
alpha(1,:) = pi0.*dens(1,:);
for t=2:T
  alpha(t,:) = (alpha(t-1,:) * A) .* dens(t,:);
  % Systematic scaling
  alpha(t,:) = alpha(t,:) /  sum(alpha(t,:));
end
if (~QUIET)
  time = cputime - time;
  fprintf(1, ' (%.2f s)\n', time);
end

% Simulation of the sequence
s = zeros(T,1);
if (~QUIET)
  fprintf(1, 'Generating random sequence...'); time = cputime;
end
% Generate last index (assumes that the last forward variable is normalized !)
s(T) = randindx(alpha(T,:),1,1);
for t=(T-1):-1:1
  % Compute probabilities conditionnal to X and s(t+1:T)
  p = A(:,s(t+1))'.*alpha(t,:);
  p = p./sum(p);
  s(t) = randindx(p,1,1);
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n',time);
end
