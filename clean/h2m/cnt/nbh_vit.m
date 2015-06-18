function [class, logl] = nbh_vit(count, TRANS, alpha, beta)

%nbh_vit   A posteriori sequence estimation for negative binomial HMM
%          using dynamic programming (also known as Viterbi algorithm).
%          Use: [class,logl] = nbh_vit(count,TRANS,alpha,beta).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 01/01/98 - 22/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Check input arguments
error(nargchk(4, 4, nargin));
% Data length
T = length(count);
if (any(count < 0) | any(count ~= fix(count)))
  error('Data does not contain positive integers.');
end
count = reshape(count, T, 1);
% Number of mixture components
N = nbh_chk(TRANS, alpha, beta);
alpha = reshape(alpha, 1, N);
beta  = reshape(beta, 1, N);

%%% 0: Compute density values on a log scale
% Compute log(count!), the second solution is usually much faster
% except if max(count) is very large
cm = max(count);
if (cm > 50000)
  dnorm = gammaln(count + 1);
else
  tmp = cumsum([0; log((1:max(count)).')]);
  dnorm = tmp(count+1);
end
logdens = ones(T,1)*(alpha.*log(beta./(1+beta)) - gammaln(alpha)) ...
    - count*log(1+beta) + gammaln(count*ones(1,N) + ones(T,1)*alpha) ...
    - dnorm*ones(1,N);

% Dynamic programming recursion
% HMM transition probabilities in log
TRANS = log(TRANS + realmin);
% Partial loglikelihood array and bactracking array
PLOGL = zeros(T, N);
BCKTR = zeros(T-1, N);
% Use uniform a priori probability for the initial state
PLOGL(1,:) = logdens(1,:) - log(N);
for t = 2:T
  [PLOGL(t,:), BCKTR(t-1,:)] = max((PLOGL(t-1,:)' * ones(1,N)) + TRANS);
   PLOGL(t,:) = PLOGL(t,:) + logdens(t,:);
end

% Backtracking
class = zeros(T,1);
[logl, class(T)] = max(PLOGL(T,:));
for t = (T-1):-1:1
  class(t) = BCKTR(t, class(t+1));
end
