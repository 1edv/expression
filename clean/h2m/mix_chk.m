function [N, p, DIAG_COV] = mix_chk (w, mu, Sigma)

%mix_chk   Checks the parameters of a mixture model and returns its dimensions.
%	Use: [N,p,DIAG_COV] = mix_chk(w,mu,Sigma).

% H2M Toolbox, Version 2.0

% Olivier Cappé, 1995 - 13/02/98
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Max deviation from 1 (beware parameters may have
% been computed in float)
MAX_DEV = 1e-5;

error(nargchk(3, 3, nargin));
% Check weight vector
[Nr, N] = size(w);
if (Nr ~= 1)
  error('Mixture weights must be specified as a row vector.');
end
if (any(w < 0) | any(w > 1))
  error('Inconsistent number in vector of mixture weights.');
end
if (abs(sum(w)-1) > MAX_DEV)
  error('Mixture weights are not normalized.');
end

% Check means of Gaussian densities
[Nr, p] = size(mu);
if (Nr ~= N)
  error('Incorrect number of mean vectors.');
end

% Check covariance matrices of Gaussian densities
[Nr, Nc] = size(Sigma);
if (Nc ~= p)
  error('The size of the covariance matrices is incorrect.');
end
if (Nr == N)
  DIAG_COV = 1;		% Test it first, for the case p = 1
elseif (Nr == p*N)
  DIAG_COV = 0;
else
  error('The size of the covariance matrices is incorrect.');
end
