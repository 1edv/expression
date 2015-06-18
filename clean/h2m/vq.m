function [CODE_n, label, dist] = vq(X, CODE, n_it, QUIET)

%vq        Vector quantization using the K-means (or LBG) algorithm.
%       Use: [CODE_n,label,dist] = vq(X,CODE,n_it)
%	Performs n_it iterations of the K-means algorithm on X, using
%	CODE as initial codebook.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 28/09/94 - 16/07/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

error(nargchk(3, 4, nargin));
if (nargin < 4)
  QUIET = 0;
end

% Dimensions of X
[n,p] = size(X);
% Codebook size
m  = length(CODE(:,1));
% Initialialize label array
label = zeros(1,n);
% As well as distortion values
dist = zeros(1,n_it);

% Main loop
CODE_n = CODE;
for iter = 1:n_it
  % 1. Find nearest neighbor for the squared distortion
  DIST = zeros(m,n);
  if (p > 1)
    for i = 1:m
      DIST(i,:) = sum(((X - ones(n,1)*CODE_n(i,:))').^2);
    end
  else
    % Beware of sum when p = 1 (!)
    DIST = (ones(m,1)*X' - CODE_n*ones(1,n)).^2;
  end
  [vm,label] = min(DIST);
  % Mean distortion
  dist(iter) = mean(vm);
  % 2. Update the codebook
  n_out = 0;
  for i = 1:m
    ind = (1:n);
    ind = ind((label == i));
    if (length(ind) == 0)
      % Isolated centroid are not modified
      n_out = n_out + 1;
    elseif (length(ind) == 1)
      % When there is only one nearest neighbor for a given codebook entry
      CODE_n(i,:) = X(ind,:);
    else
      CODE_n(i,:) = mean(X(ind,:));
    end
  end
  % Affichage
  if (~QUIET)
    fprintf(1,'Iteration %d:\t%.3f\n',iter,dist(iter));
  end
  if (n_out > 0)
    fprintf(1,'  Warning : %.0f isolated centroids\n',n_out);
  end
end
