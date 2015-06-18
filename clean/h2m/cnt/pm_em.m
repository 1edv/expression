function [wght, rate, logl,postprob] = pm_em(count, wght, rate, Nit)

%pm_em	  Estimates the parameters of a Poisson mixture using the EM algorithm.
%         Use: [wght,rate,logl,postprob] = pm_em(count,wght_0,rate_0,Nit)
%         where wght and rate are the estimated model parameters, logl
%         contains the log-likehood values for the successive
%         iterations and postprob the marginal posterior probabilities
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
N = pm_chk(wght, rate);
wght = reshape(wght, 1, N);
rate = reshape(rate, 1, N);

% Variables
logl = zeros(1,Nit);
postprob = zeros(T, N);
% Compute log(count!), the second solution is usually much faster
% except if max(count) is very large
cm = max(count);
if (cm > 50000)
  dnorm = gammaln(count + 1);
else
  tmp = cumsum([0; log((1:max(count)).')]);
  dnorm = tmp(count+1);
end

% Main loop of the EM algorithm
for nit = 1:Nit
  %%% 1: Compute a posteriori probabilities and likelihood (E)
  % Compute all densities
  postprob = exp(-ones(T,1)*rate + count*log(rate) - dnorm*ones(1,N));
  % Compute unormalized a posteriori probability
  postprob = postprob .* (ones(T,1) * wght);
  % Compute loglikelihood
  logl(nit) = sum(log(sum(postprob')));
  fprintf(1, 'Iteration %d:\t%.3f\n', (nit-1), logl(nit));
  % Normalization
  postprob = postprob ./ ((sum(postprob'))' * ones(1,N));

  %%% 2: Reestimating mixture parameters (M)
  % Unormalized weights
  wght = sum(postprob);
  % Intensities
  rate = (count.' * postprob) ./ wght;
  % Normalize weigths
  wght = wght / sum(wght);
end
