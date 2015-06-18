function [mu, Sigma, gamma] = hmm_dens (X, alpha, beta, DIAG_COV, QUIET)

%hmm_dens  Reestimates the parameters of the Gaussian distributions for an HMM.
%	Use: [mu,Sigma,gamma] = hmm_dens(X,alpha,beta,DIAG_COV).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 14/03/96 - 17/03/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Check input data
error(nargchk(4, 5, nargin));
if (nargin < 5)
  QUIET = 0;
end
if (size(alpha) ~= size(beta))
  error('Check the size of forward and backward probability matrices.');
end
[T, N] = size(alpha);
if (length(X(:,1)) ~= T)
  error('Check the number of observation vectors.');
end
p = length(X(1,:));

% Compute state occupation probabilities
if (~QUIET)
  fprintf(1, 'Reestimating gaussian parameters...'); time = cputime;
end
gamma = alpha .* beta;
gamma = gamma ./ (sum(gamma')' * ones(1,N));

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
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n', time);
end
