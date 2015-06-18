function [gamma, logl] = mix_post (X, w, mu, Sigma, QUIET)

%mix_post  A posteriori probabilities for a gaussian mixture model.
%	Use: [gamma,logl] = mix_post (X,w,mu,Sigma) returns the a posteriori
%	probabilities. logl is the log-likehood of X.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 1995 - 12/04/99
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Compiler pragmas
%#inbounds
%#realonly

% Input args.
error(nargchk(4, 5, nargin));
if (nargin < 5)
  QUIET = 0;
end
[N, p, DIAG_COV] = mix_chk(w, mu, Sigma);
[T, Nc] = size(X);
if (Nc ~= p)
  error('Observation vectors have an incorrect dimension.');
end

% Compute density values
gamma = gauseval(X, mu, Sigma, QUIET);

% A posteriori probabilities
if (~QUIET)
  fprintf(1, 'Computing a posteriori probabilities...'); time = cputime;
end
gamma = gamma .* (ones(T,1) * w);
if (nargout > 1)
  % Compute log-likelihood
  logl = sum(log(sum(gamma')));
end
gamma = gamma ./ ((sum(gamma'))' * ones(1,N));
if (~QUIET)
  time = cputime - time; fprintf(1, ' (%.2f s)\n',time);
end
