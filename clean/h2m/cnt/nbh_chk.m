function N = nbh_chk(TRANS, alpha, beta)

%nbh_chk   Checks the parameters of a negative binomial HMM
%          and returns its dimension.
%          Use: N = nbh_chk(TRANS,alpha,beta).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 31/12/97 - 17/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Constants
% Max deviation from 1 (beware parameters may have been computed in float)
MAX_DEV = 1e-6;

% Check input arguments
error(nargchk(3, 3, nargin));
% Check transition matrix
[N, Nc] = size(TRANS);
if (Nc ~= N)
  error('Transition matrix must be square.');
end
if (any(TRANS < 0) | any(TRANS > 1))
  error('Inconsistent number in transition matrix.');
end
if (any(abs(sum(TRANS')-ones(1,N)) > MAX_DEV))
  error('Transition matrix is not normalized.');
end
if (length(alpha) ~= N)
  error('Transition matrix and shape parameters must have the compatible sizes.');
end
if (any(alpha) <= 0)
  error('Shape parameters must be positive.');
end
if (length(beta) ~= N)
  error('Transition matrix and inverse scales must have the compatible sizes.');
end
if (any(beta) <= 0)
  error('Inverse scales must be positive.');
end
