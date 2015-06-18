function [TRANS, rate, logl, beta, dens] = ph_em(count, TRANS, rate, Nit)

%ph_em	  Estimates the parameters of a Poisson HMM using the EM algorithm.
%         Use: [TRANS,rate,logl,postprob,dens] =
%         ph_em(count,TRANS_0,rate_0,Nit). where TRANS and rate are
%         the estimated model parameters, logl contains the
%         log-likehood values for the successive iterations, postprob
%         the marginal posterior probabilities at the last iteration
%         and dens the negative binomial probabilities computed also
%         at the last iteration.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 30/12/97 - 16/08/2001
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
N = ph_chk(TRANS, rate);
rate = reshape(rate, 1, N);

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
alpha = zeros(T, N);
scale = [1 ; zeros(T-1, 1)];
beta = zeros(T, N);

% Main loop of the EM algorithm
for nit = 1:Nit
  %%% 1: E-Step, compute density values
  dens = exp(-ones(T,1)*rate + count*log(rate) - dnorm*ones(1,N));
  
  %%% 2: E-Step, forward recursion and likelihood computation
  % Use a uniform a priori probability for the initial state
  alpha(1,:) = dens(1,:)/N;
  for t = 2:T
    alpha(t,:) = (alpha(t-1,:) * TRANS) .* dens(t,:);
    % Systematic scaling
    scale(t) = sum(alpha(t,:));
    alpha(t,:) = alpha(t,:) / scale(t);
  end
  % Compute log-likelihood
  logl(nit) = sum(log(scale));
  fprintf(1, 'Iteration %d:\t%.3f\n', (nit-1), logl(nit));

  %%% 3: E-Step, backward recursion
  % Scale the backward variable with the forward scale factors (this ensures
  % that the reestimation of the transition matrix below is correct)
  beta(T, :) = ones(1, N);
  for t = (T-1):-1:1
    beta(t, :) = (beta(t+1,:).* dens(t+1,:)) * TRANS';
    % Apply scaling
    beta(t,:) = beta(t,:) / scale(t);
  end
  
  %%% 4: M-Step, reestimation of the transition matrix
  % Compute unnormalized transition probabilities (this is indeed still the
  % end of the E-step, which explains that TRANS appears on the right-hand
  % side below)
  TRANS = TRANS .* (alpha(1:(T-1),:).' * (dens(2:T,:) .* beta(2:T,:)));
  % Normalization of the transition matrix
  TRANS = TRANS ./ (sum(TRANS.').' *ones(1,N));
  
  %%% 5: M-Step, reestimation of the rates
  % Compute a posteriori probabilities (and store them in matrix beta
  % to save some space)
  beta = alpha .* beta;
  beta = beta ./ (sum(beta.').' * ones(1,N));
  % Reestimate rates
  rate = (count.' * beta) ./ sum(beta);
end
