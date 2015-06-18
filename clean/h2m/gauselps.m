function gauselps(mu, Sigma, lt, level, c)

%gauselps  Plots 2D projection of gaussian ellipsoids.
%	Use: gauselps(mu,Sigma[,lt,level,c]) . Optional arguments:
%	lt: line type, level: probability inside the ellipse (default is
%	0.39 for standard deviation, other values include 0.5, 0.75, 0.9,
%	0.99), c = [c1 c2]: coordinates to be considered (default is [1 2]).

% H2M Toolbox, Version 2.0
% Olivier Cappé, 1994 - 12/05/99
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Input arguments
error(nargchk(2, 5, nargin));
if (nargin < 5)
  c = [1 2];
  if (nargin < 4)
    % Plot the standard deviation on each axis by default
    level = 0.39;
    if (nargin < 3)
      lt = '-';
    end
  end
end
[N,p] = size(mu);
[nc, nr] = size(Sigma);
if (nr ~= p)
  error('Dimension of the covariance matrices is incorrect.');
end
if (nc == N)
  DIAG_COV = 1;
elseif (nc == N*p)
  DIAG_COV = 0;
else
  error('Dimension of the covariance matrices is incorrect.');
end
if (p == 1)
  error('Can''t plot ellipsoids for scalar mixture!');
end
% Level of the ellipse
tab_chi2 = [0.39 1; 0.5 1.177; 0.75 1.665; 0.9 2.146; 0.95 2.448; 0.99 3.035];
% Result will be extremely poor if level is not in the table...
hom = table1(tab_chi2, level);

% 90% circle
nt = 50;
theta = linspace(0, 2*pi-eps, nt);
circle = hom * [cos(theta) ; sin(theta)];

% Plot
hold_state = ishold;
if (DIAG_COV)
  for i=1:N
    elps = mu(i,c)' * ones(1,nt) + diag(sqrt(Sigma(i,c))) * circle;
    plot(elps(1,:), elps(2,:), lt);
    hold on;
    text(mu(i,c(1)), mu(i,c(2)), int2str(i));
  end
else
  for i =1:N
    F = chol(Sigma(c+(i-1)*p,c))';
    elps = mu(i,c)' * ones(1,nt) + F * circle;
    plot(elps(1,:), elps(2,:), lt);
    hold on;
    text(mu(i,c(1)), mu(i,c(2)), int2str(i));
  end
end
% Restore hold mode
if (~hold_state)
  hold off;
end
