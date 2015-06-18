function [alpha, beta] = nb_ml(count, TOL)

%nb_ml	  Maximum likelihood estimates for negative binomial data.
%         Use: [alpha, beta] = nb_ml(count[, TOL]).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 05/02/98 - 22/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Known problems: The optimization fails when there are very few data points
% and/or when most of them are zero (one example of problematic data set is
% count = [0 1 0 0])

% Constants
BETA_INIT_DEF = 10;
NIT_MAX = 250;

% Check input arguments
error(nargchk(1, 2, nargin));
if (nargin < 2)
  % Tolerance on the gradient for stopping the iterations
  TOL = 1e-5;
end
T = length(count);

% Start with moment estimates
m1 = mean(count);
beta = 1./(cov(count)/m1-1);
if (beta < 0)
  fprintf(1, 'Warning: N-Binomial model doesn''t seem to fit the data\n');
  beta = BETA_INIT_DEF;
end
alpha = m1*beta;

% Use coordinate ascent method (which is globaly convergent here, because
% beta/alpha has a unique maximizer and alpha/beta is convex)
STOP = 0;
it_count = 0;
while (~STOP)
  % Update beta
  beta = alpha/m1;

  % Use a Newton step for updating alpha
  % Compute first derivative
  grad = T*(log(beta/(1+beta)) - digamma(alpha)) + sum(digamma(count+alpha));
  % and second derivative
  hess = -T*trigamma(alpha) + sum(trigamma(count+alpha));
  % Newton step
  step = grad/hess;
  alpha = alpha - step;

  % When performing the Newton step, one should check that the likelihood
  % indeed increases and that alpha does not become negative (in practise this
  % is never needed - except in the extreme cases mentioned above).

  % Check stopping condition on alpha only
  STOP = ((abs(step/alpha) < TOL) | (it_count >= NIT_MAX)) ;
  it_count = it_count + 1;
end
if (it_count > NIT_MAX)
  fprintf(1, 'Warning: Maximum number of iterations (%d) exceeded\n', NIT_MAX);
else  
  fprintf(1, '\n%d Iterations.\n', it_count);
end
