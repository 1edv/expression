%ex_bic    Example of cluster analysis with varying number of mixture
%          components using the BIC (Bayesian Information Criterion) to assess
%          the fit provided by mixture model of different sizes.
%
%       The funny thing here is to play with the parameters of the simulated
%       data (mixture parameters and number of outliers) and with the
%       estimation parameters (initialization procedure and number of
%       iterations). It turns out that the decision is somewhat sensitive to
%       these (although it usually gives 3)...
%
%	Beware : This is not a real demo, so you should take a look at the
%       source file before executing it if you wan't to understand what's
%       going on! If using octave --traditional, you should also add the
%       subdirectory h2m/octave to your loadpath using the path command.

% H2M Toolbox, Version 2.0
% O. Cappe, 19/05/99 - 24/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

clear;

% Constant: Flag for full covariance matrices
FULL_COV = 0;

% Generate some 2D synthetic data (Gaussian mixture model with 3 compoenent
% and some outliers)
% Total number of observations
T = 200;
% Number of outliers
OUTL_NB = 5;
% Mixture data
% Data dimension
p = 2;
% Number of components
N = 3;
w = [0.4 0.4 0.2];
mu = [1 0; 0 1; -1 0];
Sigma = [0.1 0; 0 0.1; 0.1 -0.03; -0.03 0.1; 0.2 0; 0 0.2];
X = mix_gen(w, mu, Sigma, T-OUTL_NB);
% Append outliers (this of course doesn't look like iid data anymore but the
% following analysis would be exactly the same if the data points had been
% shuffled randomly)
X = [X; 4*(rand(OUTL_NB,2)-0.5)];

% Cluster analysis
% Maximal size of mixture model
NC_MAX = 5;
% Values of the BIC criterion
BIC = zeros(1,NC_MAX);
% Number of EM iterations
NIT = 20;

% NC = 1 corresponds to a gaussian model
NC = 1;
w = 1;
mu = mean(X);
Sigma = cov(X);
logl = sum(gauslogv(X, mu, Sigma,1));
BIC(NC) = 2*logl-(p*NC+p*NC+NC-1)*log(T);

% Plot data fit
clf;
% Data points
subplot(3, 3, 1);
plot(X(:,1), X(:,2), '.');
axis([-2 2 -2 2]);
% Gausian elipsoid with 50% of the data
subplot(3, 3, 2);
title(['BIC = ' sprintf('%.0f', BIC(NC))]);
gauselps(mu,Sigma,'-',0.5);
axis([-2 2 -2 2]);
drawnow;

for NC = 2:NC_MAX
  NC
  % Here is an example of a simple heuristic for mixture-based hierarchic
  % clustering which may be viewed as an alternative to vector quantization for
  % the initialization of the model (it also has the advantage that one is not
  % restricted to powers of 2 for NC) Find most represented class
  [tmp, ind] = max(w);
  % Find main axis
  [V, D] = eig(Sigma(1+(ind-1)*p:ind*p,:));
  [tmp, ind2] = max(diag(D));
  % Duplicate one component
  mu = [mu(1:ind-1,:); mu(ind,:)-sqrt(D(ind2,ind2))*V(:,ind2)'; ...
	mu(ind,:)+sqrt(D(ind2,ind2))*V(:,ind2)'; mu(ind+1:NC-1,:)];
  Sigma = [Sigma(1:ind*p,:); Sigma(1+(ind-1)*p:(NC-1)*p,:)];
  % EM iterations
  w = [w(1:ind-1) w(ind)/2 w(ind)/2 w(ind+1:NC-1)];
  % Compute log-likelihood
  [w, mu, Sigma, logl] = mix(X, w, mu, Sigma, NIT);
  [gamma, logl] = mix_post(X, w, mu, Sigma, 1);
  % and BIC criterion
  BIC(NC) = 2*logl-(p*NC+p*NC+NC-1)*log(T);
  % Plot data fit
  subplot(3, 3, 1+NC);
  % Gausian elipsoids with 50% of the data
  title(['BIC = ' sprintf('%.0f', BIC(NC))]);
  gauselps(mu,Sigma,'-',0.5);
  axis([-2 2 -2 2]);
  drawnow;
end

% When evrything is done report BIC criterion (best fit is given by the
% maximal value)
subplot(3,1,3);
axis;
plot(BIC);
title('BIC');
BIC
