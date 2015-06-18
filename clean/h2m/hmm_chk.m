function [N, p, DIAG_COV] = hmm_chk (A, pi0, mu, Sigma)

%hmm_chk   Checks the parameters of an HMM and returns its dimensions. 
%	Use : [N,p,DIAG_COV] = hmm_chk(A,pi0,mu,Sigma).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 1995 - 13/02/98
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Max deviation from 1 (beware parameters may have
% been computed in float)
MAX_DEV = 1e-5;

error(nargchk(4, 4, nargin));
% Check transition matrix
[N, Nc] = size(A);
if (Nc ~= N)
  error('Transition matrix must be square.');
end
if (any(A < 0) | any(A > 1))
  error('Inconsistent number in transition matrix.');
end
if (any(abs(sum(A')-ones(1,N)) > MAX_DEV))
  error('Transition matrix is not normalized.');
end

% Check initial probability vector
[Nr, Nc] = size(pi0);
if (Nr ~= 1)
  error('Initial probabilities must be specified as a row vector.');
end
if (Nc ~= N)
  error('Vector of initial probabilities has an incorrect size.');
end
if (any(pi0 < 0) | any(pi0 > 1))
  error('Inconsistent number in initial probabilities vector.');
end
if (abs(sum(pi0)-1) > MAX_DEV)
  error('Initial probabilities are not normalized.');
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
  DIAG_COV = 1;		% Test it first for the case p = 1
elseif (Nr == p*N)
  DIAG_COV = 0;
else
  error('The size of the covariance matrices is incorrect.');
end
