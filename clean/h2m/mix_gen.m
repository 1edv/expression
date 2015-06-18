function [X, s] = mix_gen (w, mu, Sigma, T)

%mix_gen   Generates a sequence of observation for a gaussian mixture model.
%	Use: [X,s] = mix_gen(w,mu,Sigma,T).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 1995 - 05/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (~((exist('randindx') == 2) | (exist('randindx') == 3)))
  error('Function randindx is missing.');
end

error(nargchk(4, 4, nargin));
[N, p, DIAG_COV] = mix_chk(w, mu, Sigma);

% Generate sequence of indexes
s = randindx(w, T);

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
