function [CODE, label, dist] = svq(X, lev, n_it)

%svq	  Vector quantization using successive binary splitting steps.
%	Use: [CODE,label,dist] = svq(X,lev,n_it).
%	The final codebook dimension lev should be a power of two. dist
%	returns the distorsion values at the end of intermediate step.
%	n_it is the number of iterations performed in each step.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 28/09/94 - 04/03/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (exist('vq') ~= 2)
  error('Function vq is missing.');
end
% Turn verbose mode off
QUIET = 1;
% Input agruments
error(nargchk(3, 3, nargin));
% Dimension of imput data
[n,p] = size(X);
% Number of spliting steps
nbs = round(log2(lev));
lev = 2^nbs;
% Fixed perturbation
perturb = 0.01;

% Initialize first centroid with global mean
CODE = zeros(lev, p);
CODE_ = zeros(lev, p);
CODE(1,:) = mean(X);
label = ones(n,1);

for i=1:nbs
  % 1. Codebook splitting
  for j=1:(2^(i-1))
    CODE_(2*j-1,:) = (1+perturb) * CODE(j,:);
    CODE_(2*j,:)   = (1-perturb) * CODE(j,:);
  end
  % 2. K-means optimization
  [CODE(1:2^i,:),label,vdist] = vq(X,CODE_(1:2^i,:),n_it,QUIET);
  dist(i) = vdist(n_it);
  fprintf(1, 'Codebook size %d:\t%.3f\n',2^i,dist(i));
end
