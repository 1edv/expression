%ex_basic  How to use the h2m functions on the three basic model types:
%          1. Fully connected HMM 2. Left-right HMM 3. Mixture model
%
%	Beware : This is not a real demo, so you should take a look at the
%       source file before executing it if you wan't to understand what's
%       going on! If using octave --traditional, you should also add the
%       subdirectory h2m/octave to your loadpath using the path command.

% H2M Toolbox, Version 2.0
% O. Cappe, 24/06/96 - 14/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% The notations are essentially those of
%   AUTHOR   = "L. R. Rabiner",
%   TITLE    = "A tutorial on hidden Markov models and selected applications
%              in speech recognition",
%   JOURNAL  = procieee,
%   VOLUME   = 77,
%   NUMBER   = 2,
%   MONTH    = feb,
%   PAGES    = "257-285",
%   YEAR     = 1989
%
% pi0:	 Vector of initial probabilities, pi_0(i) = P(s_1 = i)
% A:	 Transition matrix a(i,j) = P(s_{t+1}=j | s_{t}=i)
% mu:    Mean vector of gaussian densities (stacked one above the other)
% Sigma: Covariance matrices (or diagonal of) stacked one above the other
%
% See the on-line help of the functions for more information.

clear;

% Constant: Flags for diagonal and full covariance matrices
FULL_COV = 0;
DIAG_COV = 1;

% Example 1: Ergodic model, with full covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data generation
p = 2;
N = 3;
T = 100;
n_it = 10;
pi0 = rand(1,N); pi0 = pi0/ sum(pi0);
A = rand(N); A = A ./ (sum(A')' * ones(1,N));
mu = [1 0; 0 1; -1 0];
Sigma = [0.1 0; 0 0.1; 0.1 -0.03; -0.03 0.1; 0.2 0; 0 0.2];
X = hmm_gen(A, pi0, mu, Sigma, T);
% Initialization values
pi0_ = pi0;
A_ = A;
mu_ = [0.2 0; 0 0; -0.2 0];
Sigma_ = 0.2*[diag(ones(1,p)); diag(ones(1,p)); diag(ones(1,p))];
% EM iterations
for i = 1:n_it
  [alpha, beta, logl_tmp, dens] = hmm_fb(X, A_, pi0_, mu_, Sigma_);
  logl(i) = logl_tmp;
  [A_, pi0_] = hmm_tran(alpha, beta, dens, A_, pi0_);
  [mu_, Sigma_] = hmm_dens(X, alpha, beta, FULL_COV);
end
% You could use more simply
% [A_,pi0_,mu_,Sigma_,logl] = hmm(X,A_,pi0_,mu_,Sigma_,n_it);
% Plot estimated means and covariances
clf;
gauselps(mu_,Sigma_);
hold on;
gauselps(mu,Sigma, 'r');

pause;

% Example 2: Left-right model, with diagonal covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data generation
p = 2;
N = 3;
T = 30;
K = 10;
n_it = 10;
pi0 = sparse([1 0 0]);
A = sparse([0.8 0.16 0.04; 0 0.9 0.1; 0 0 1]);
mu = [1 0; 0 1; -1 0];
Sigma = [0.1 0.1; 0.1 0.2; 0.2 0.2];
for i = 1:K
  X((1+(i-1)*T:i*T),:) = hmm_gen(A, pi0, mu, Sigma, T);
end
st = 1 + (0:(K-1))*T;
% Initialization values
A_ = sparse([0.4 0.4 0.2; 0 0.5 0.5; 0 0 1]);;
[mu_, Sigma_] = hmm_mint(X, st, N, DIAG_COV);
% EM iterations
for i = 1:n_it
  [A_, logl(i), gamma] = hmm_mest(X, st, A_, mu_, Sigma_);
  [mu_,Sigma_] = mix_par(X, gamma, DIAG_COV);
end
% You could use more simply
% [A_,mu_,Sigma_,logl] = lrhmm (X,st,A,mu,Sigma,n_it);
% Plot estimated means and covariances
clf;
gauselps(mu_,Sigma_);
hold on;
gauselps(mu,Sigma, 'r');

pause;

% Example 3: Gaussian Mixture Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data generation
p = 2;
N = 4;
T = 500;
n_it = 10;
w = [0.2 0.4 0.2 0.2];
mu = [1 0; 0 1; -1 0; 0 -1];
Sigma = [0.1 0; 0 0.1; 0.1 -0.03; -0.03 0.1; 0.2 0; 0 0.2; 0.1 0; 0 0.1];
X = mix_gen(w, mu, Sigma, T);
% Initial values (using svq)
w_ = ones(1,N)/N;
[mu_, label] = svq(X,N,10);
Sigma_ = [cov(X((label == 1),:)); cov(X((label == 2),:)); cov(X((label == 3),:)); cov(X((label == 4),:))];
% EM iterations
for i = 1:n_it
  [gamma, logl(i)] = mix_post(X, w_, mu_, Sigma_);
  [mu_, Sigma_, w_] = mix_par(X, gamma, FULL_COV);
end
% You could use more simply
% [w_,mu_,Sigma_,logl] = mix (X,w,mu,Sigma,n_it);
% Plot estimated means and covariances
clf;
gauselps(mu_, Sigma_);
hold on;
gauselps(mu, Sigma, 'r');
