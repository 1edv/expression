function [TRANS, alpha, beta, logl, bckwrd,dens] = nbh_em(count, TRANS, alpha, beta, Nit)

%nbh_em	  Estimates the parameters of a negative binomial HMM using EM.
%         Use: [TRANS,alpha,beta,logl,postprob,dens] =
%         nbh_em(count,TRANS_0,alpha_0,beta_0,Nit) where TRANS, alpha
%         and beta are the estimated model parameters, logl contains
%         the log-likehood values for the successive iterations,
%         postprob the marginal posterior probabilities at the last
%         iteration and dens the negative binomial probabilities
%         computed also at the last iteration.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 31/12/97 - 22/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Note: This actually uses the ECM (or SAGE) version of EM since
% complete maximization of the EM intermediate quantity is not
% directly feasible, but maximization with respect to the beta and
% alpha coordinates is possible (using Newton algorithm for the
% latter).
%
% Refs.:
% - X. L. Meng, D. B. Rubin, Maximum likelihood estimation via the
%   ECM algorithm: A general framework, Biometrika, 80(2):267-278 (1993).
% - J. A. Fessler, A. O. Hero, Space-alternating generalized
%   expectation-maximization algorithm, IEEE Tr. on Signal
%   Processing, 42(10):2664 -2677 (1994).

% Known problems: The optimization fails when there are very few data points
% and/or when most of them are zero. Using extremely bad starting values for
% alpha and beta can also cause the optimization to fail. In such cases, the
% likelihood may even decrease from one iteration to the other.

% Check input arguments
error(nargchk(5, 5, nargin));
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

% Compute log(count!), the second solution is usually much faster
% except if max(count) is very large
cm = max(count);
if (cm > 50000)
  dnorm = gammaln(count + 1);
else
  tmp = cumsum([0; log((1:max(count)).')]);
  dnorm = tmp(count+1);
end

% Variables
logl = zeros(1,Nit);
dens = zeros(T, N);
scale = [1 ; zeros(T-1, 1)];
% The forward and backward variables
forwrd = zeros(T, N);
bckwrd = zeros(T, N);

% Main loop of the EM algorithm
for nit = 1:Nit
  %%% 1: E-Step, compute density values
  dens = exp(ones(T,1)*(alpha.*log(beta./(1+beta)) - gammaln(alpha)) ...
      - count*log(1+beta) + gammaln(count*ones(1,N) + ones(T,1)*alpha) ...
      - dnorm*ones(1,N));
  
  %%% 2: E-Step, forward recursion and likelihood computation
  % Use a uniform a priori probability for the initial state
  forwrd(1,:) = dens(1,:)/N;
  for t = 2:T
    forwrd(t,:) = (forwrd(t-1,:) * TRANS) .* dens(t,:);
    % Systematic scaling
    scale(t) = sum(forwrd(t,:));
    forwrd(t,:) = forwrd(t,:) / scale(t);
  end
  % Compute log-likelihood
  logl(nit) = log(sum(forwrd(T,:))) + sum(log(scale));
  fprintf(1, 'Iteration %d:\t%.3f\n', (nit-1), logl(nit));

  %%% 3: E-Step, backward recursion
  % Scale the backward variable with the forward scale factors (this ensures
  % that the reestimation of the transition matrix is correct)
  bckwrd(T, :) = ones(1, N);
  for t = (T-1):-1:1
    bckwrd(t, :) = (bckwrd(t+1,:).* dens(t+1,:)) * TRANS';
    % Apply scaling
    bckwrd(t,:) = bckwrd(t,:) / scale(t);
  end
  
  %%% 4: M-Step, reestimation of the transition matrix
  % Compute unnormalized transition probabilities (this is indeed still the
  % end of the E-step, which explains that TRANS appears on the right-hand
  % side below)
  TRANS = TRANS .* (forwrd(1:(T-1),:).' * (dens(2:T,:) .* bckwrd(2:T,:)));
  % Normalization of the transition matrix
  TRANS = TRANS ./ (sum(TRANS.').' *ones(1,N));
  
  %%% 5: CM-Step 1, reestimation of the inverse scales beta with alpha fixed
  % Compute a posteriori probabilities (and store them in matrix
  % bckwrd to save some space)
  bckwrd = forwrd .* bckwrd;
  bckwrd = bckwrd ./ (sum(bckwrd.').' * ones(1,N));
  % Reestimate shape parameters beta conditioning on alpha
  eq_count = sum(bckwrd);
  beta = alpha./((count.'*bckwrd)./eq_count);

  %%% 5: CM-Step 2, reestimation of the shape parameters with beta fixed
  % Use digamma and trigamma function to perfom a Newton step on
  % the part of the intermediate quantity of EM that depends on alpha
  % Compute first derivative for all components
  grad = eq_count .* (log(beta./(1+beta)) - digamma(alpha)) ...
    + sum(bckwrd .* digamma(count*ones(1,N) + ones(T,1)*alpha));
  % And second derivative
  hess = -eq_count.*trigamma(alpha) ...
    + sum(bckwrd .* trigamma(count*ones(1,N) + ones(T,1)*alpha));
  % Newton step
  tmp_step = - grad./hess;
  tmp = alpha + tmp_step;
  
  % When performing the Newton step, one should check that the intermediate
  % quantity of EM indeed increases and that alpha does not become negative. In
  % practise this is almost never needed but the code below may help in some
  % cases (when using real bad initialization values for the parameters for
  % instance)
  while (any(tmp <= 0))
    fprintf(1, 'Warning: could not update alpha\n');
    tmp_step = tmp_step/10;
    tmp = alpha + tmp_step;
  end
  alpha = tmp;
end
