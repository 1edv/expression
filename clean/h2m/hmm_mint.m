function [mu, Sigma] = hmm_mint (X, st, N, DIAG_COV, QUIET)

%hmm_mint  Initializes the model parameters using multiple observations
%          for a left-to-right HMM.
%	Use: [mu,Sigma] = hmm_mint(X,st,N,DIAG_COV).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 1995 - 25/11/99
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Check input data
error(nargchk(4, 5, nargin));
if (nargin < 5)
  QUIET = 0;
end

% Dimensions
[KT, p] = size(X);
K = length(st);
% Check the data sections info
if ((st(1) ~= 1) | (st(K) >= KT))
  error('Check the limits of the observation sequences.');
end
% Check that there is enough data: The criterion below is based on the fact
% that in the diagonal case, fewer than 2 observations results in a null
% estimate for the variances, and in the full covariance case, less than p
% observations implies a singular covariance matrix estimate (this is certainly
% not very conservative and does not guarantee in any way that the estimate
% will be "correct").
if (DIAG_COV)
  if (any((diff(st)*K/N) < 2))
    error('Too many states or not enough data (singular variance estimation).');
  end
else
  if (any((diff(st)*K/N) < p))
    error('Too many states or not enough data (singular covariance estimation).');
  end
end
% In any case one shouldn't accept segments shorther than 2 vectors
if (any(round((diff(st)/N)) < 2))
  error('Too many states or not enough data (one of the segments contains less than 2 observations).');
end

if (~QUIET)
  fprintf(1, 'Initializing gaussian parameters...'); time = cputime;
end
mu = zeros(N,p);
count = zeros(1, N);
if (DIAG_COV)
  % Initialization
  Sigma = zeros(N,p);
  for i=1:K
    % Limits of the observation sequence
    first = st(i);
    if (i == K)
      last = KT;
    else
      last = st(i+1)-1;
    end
    % Uniform segmentation
    seg = round(linspace(first, last, N+1));
    % Frame count
    count = count + (seg(2:N+1)-seg(1:N)+1);
    % Estimate for each state
    for j=1:N
      mu(j,:) = mu(j,:) + sum(X(seg(j):seg(j+1),:));
      Sigma(j,:) = Sigma(j,:) + sum(X(seg(j):seg(j+1),:).^2);
    end
  end
  % Normalization
  mu = mu ./ (count' * ones(1,p));
  Sigma = Sigma ./ (count' * ones(1,p)) - mu.^2;
else
  % Initialization
  Sigma = zeros(N*p,p);
  for i=1:K
    % Limits of the observation sequence
    first = st(i);
    if (i == K)
      last = KT;
    else
      last = st(i+1)-1;
    end
    % Uniform segmentation
    seg = round(linspace(first, last, N+1));
    % Frame count
    count = count + (seg(2:N+1)-seg(1:N)+1);
    % Estimate for each state
    for j=1:N
      mu(j,:) = mu(j,:) + sum(X(seg(j):seg(j+1),:));
      Sigma(1+(j-1)*p:j*p,:) = Sigma(1+(j-1)*p:j*p,:)...
        + X(seg(j):seg(j+1),:)'*X(seg(j):seg(j+1),:);
    end
  end
  % Normalization
  mu = mu ./ (count' * ones(1,p));whos
  for j = 1:N
    Sigma(1+(j-1)*p:j*p,:) = Sigma(1+(j-1)*p:j*p,:) / count(j)...
      - mu(j,:)' * mu(j,:);
  end
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n', time);
end
