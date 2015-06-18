function [X, s] = hmm_gen (A, pi0, mu, Sigma, T)

%hmm_gen   Generates a sequence of observation given a HMM
%	   with gaussian state-conditional densities.
%	Use : [X,s] = hmm_gen(A,pi0,mu,Sigma,T|s). If the
%	last argument is of length one, it is taken as the
%	number of observations, otherwise it is considered
%	as a specified state sequence (which means that s=T).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 11/03/96 - 05/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (~((exist('randindx') == 2) | (exist('randindx') == 3)))
  error('Function randindx is missing.');
end

% Imput arguments
error(nargchk(5, 5, nargin));
[N, p, DIAG_COV] = hmm_chk(A, pi0, mu, Sigma);

% State sequence
if (length(T) == 1)
  s = zeros(T,1);
  % Generate state sequence
  for i = 1:T
    % Generate state sequence
    if (i == 1)
      s(i) = randindx(pi0, 1, 1);
    else
      s(i) = randindx(A(s(i-1), :), 1, 1);
    end
  end
else
  s = T;
  T = length(s);
  % Minimal check
  if (any(s > N) | any(s < 1))
    error('Last argument is not a valid state sequence');
  end
end

% Generate random numbers all at once
X = randn(T, p);
if (DIAG_COV)
  % Generate observations
  X = mu(s, :) + X .* sqrt(Sigma(s, :));
else
  % Cholevsky decomposition for full matrices
  F = zeros(size(Sigma));
  for i=1:N
    % see use of F below
    F((1+(i-1)*p):(i*p),:) = chol(Sigma((1+(i-1)*p):(i*p),:));
  end
  % Main loop
  for i = 1:T
    % Generate observation vector
    % This is OK since Sigma = F'*F
    X(i, :) =  mu(s(i), :) + X(i,:) * F((1+(s(i)-1)*p):(s(i)*p),:);
  end
end
