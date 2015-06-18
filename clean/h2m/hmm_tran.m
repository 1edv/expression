function [A, pi0] = hmm_tran (alpha, beta, dens, A_prev, pi0_prev, QUIET)

%hmm_tran  Reestimation the transition matrix and initial dist. of an HMM.
%	Use : [A,pi0] = hmm_tran(alpha,beta,dens,A_prev,pi0_prev). If A_prev
%	and piO_prev are sparse, the output estimates are also sparse (and
%	only the relevant values are reestimated).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 12/03/96 - 30/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Dimensions (without checking A_prev and pi0_prev)
error(nargchk(5, 6, nargin));
if (nargin < 6)
  QUIET = 0;
end
if ((size(alpha) ~= size(beta)) | (size(alpha) ~= size(dens)))
  error('First three matrices should have the same size.');
end
[T, N] = size(alpha);

% Initial probabilities
if (~QUIET)
  fprintf(1, 'Reestimating transition parameters...'); time = cputime;
end
if (issparse(pi0_prev))
  ind = find(pi0_prev);
  pi0 = sparse(zeros(1,N));
  pi0(ind) = alpha(1,ind) .* beta(1,ind) / sum(alpha(1,ind) .* beta(1,ind));
else
  pi0 = alpha(1,:) .* beta(1,:) / sum(alpha(1,:) .* beta(1,:));
end

% Transition matrix
if (issparse(A_prev))
  A = sparse(zeros(size(A_prev)));
  for i=1:N
    ind = find(A_prev(i,:));
    num = (alpha(1:(T-1),i) * ones(1,length(ind))) .* dens(2:T,ind) .* beta(2:T,ind);
    A(i,ind) = A_prev(i,ind) .* sum(num);
    A(i,ind) = A(i,ind) / sum(A(i,ind));
  end
else
  A = A_prev .* (alpha(1:(T-1),:).' * (dens(2:T,:) .* beta(2:T,:)));
  % Normalization of the transition matrix
  A = A ./ (sum(A.').' *ones(1,N));
end
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n', time);
end
