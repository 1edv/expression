function w = statdis(A)

%statdis   Returns the stationary distribution of a Markov chain.
%	Use : w = statdis(A), where A is the transition matrix.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 16/07/97 - 18/07/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Max deviation from 1 (beware parameters may have
% been computed in float)
MAX_DEV = 1e-6;

% Input args
error(nargchk(1, 1, nargin));

% Check that A is a transition matrix
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

% Compute left eigenvalues
[V,D] = eig(A');
d = diag(D);
n1 = (d > 1-MAX_DEV);
if (sum(n1) == 1)
  % Stationnary distribution is given by the left eigenvector corresponding
  % to the eigenvalue 1
  w = V(:,n1)';
  w = w./sum(w);
elseif (sum(n1) > 1)
  fprintf(1, 'Warning: transition matrix is not irreducible\n');
  w = V(:,n1)';
  for i = 1:length(w(:,1))
    w(i,:) =  w(i,:)./sum(w(i,:));
  end
else
  % This should never happen if A is a transition matrix
  error('matrix does not seem to be a transition matrix');
end
