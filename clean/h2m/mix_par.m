function [mu, Sigma, w] = mix_par (X, gamma, DIAG_COV, QUIET)

%mix_par   Reestimate mixture parameters.
%	Use: [mu,Sigma,w] = mix_par(X,gamma,DIAG_COV).
%	Note that mix_par can also be used for re-estimating HMM parameters
%	from posterior state probabilities (w is omitted in this case).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 25/03/96 - 12/04/99
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Compiler pragmas
%#inbounds
%#realonly

% Check input data
error(nargchk(3, 4, nargin));
if (nargin < 4)
  QUIET = 0;
end
[T, N] = size(gamma);
if (length(X(:,1)) ~= T)
  error('Check the number of observation vectors.');
end
p = length(X(1,:));

if (~QUIET)
  fprintf(1, 'Reestimating mixture parameters...'); time = cputime;
end
mu = zeros(N, p);
if (DIAG_COV)
  Sigma = zeros(N, p);
  for i=1:N
    mu(i,:) = sum (X .* (gamma(:,i) * ones(1,p)));
    Sigma(i,:) = sum(X.^2 .* (gamma(:,i) * ones(1,p)));
  end
  % Normalization
  mu = mu ./ (sum(gamma)' * ones(1,p));
  Sigma = Sigma ./ (sum(gamma)' * ones(1,p));
  % Sigma
  Sigma = Sigma - mu.^2;
else
  Sigma = zeros(N*p, p);
  for i = 1:N
    for j = 1:T
      mu(i,:) = mu(i,:) + gamma(j,i)*X(j,:);
      Sigma(((1+(i-1)*p):(i*p)),:) = Sigma(((1+(i-1)*p):(i*p)),:)...
        + gamma(j,i)*(X(j,:)'*X(j,:));
    end
    % Normalization
    mu(i,:) = mu(i,:)/sum(gamma(:,i));
    Sigma(((1+(i-1)*p):(i*p)),:) = Sigma(((1+(i-1)*p):(i*p)),:) /sum(gamma(:,i));
    % Sigma
    Sigma(((1+(i-1)*p):(i*p)),:) = Sigma(((1+(i-1)*p):(i*p)),:)...
        - mu(i,:)'*mu(i,:);
  end
end
% Mixture proportions
if (nargout > 2)
  w = sum(gamma);
  w = w / sum(w);
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n', time);
end
