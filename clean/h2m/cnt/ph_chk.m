function N = ph_chk(TRANS, rate)

%ph_chk	  Checks the parameters of a Poisson HMM and returns its dimension.
%         Use: N = ph_chk(TRANS,rate).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 30/12/97 - 31/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Constants
% Max deviation from 1 (beware parameters may have been computed in float)
MAX_DEV = 1e-6;

% Check input arguments
error(nargchk(2, 2, nargin));
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
if (length(rate) ~= N)
  error('Transition matrix and rates must have the compatible sizes.');
end
if (any(rate) <= 0)
  error('Poisson rates must be positive.');
end
