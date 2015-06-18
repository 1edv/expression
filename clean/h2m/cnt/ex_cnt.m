%ex_cnt    Script to illustrate the three basic models handled by H2m/cnt
%            - Poisson mixtures
%            - Poisson hidden Markov models
%            - Negative-binomial hidden Markov models
%          using synthetic data.
%
%	Beware : This is not a real demo, so you should take a look at the
%       source file before executing it if you wan't to understand what's
%       going on!

% H2M/cnt Toolbox, Version 2.0
% Olivier Capp�, 30/12/97 - 17/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

clear all;
path(path, '..');

% Try to figure out which of matlab or octave is used (this is crude)
s_ver = version;
if (s_ver(1) == '2')
  S_VER = 0; % Octave >= 2.0.14
  fprintf(1, ['I am assuming that your are using GNU Octave version\n' ...
    '2.0.14 or above in --traditional mode (congratulations!).\n']);
else
  S_VER = 1;
  fprintf(1, ['I am assuming that your are using Matlab with\n' ...
    'the Statistics toolbox.']);
end
tmp = input(['If you are not happy with that enter one of\n  [0] GNU' ...
  ' Octave --traditional\n  [1] Matlab\nor just' ...
  ' hit <return> if you believe that the current choice is right : ']);
if (~isempty(tmp))
  S_VER = tmp;
end
if (S_VER == 0)
  path(path, '../octave');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '\nPoisson mixture model.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
T = 80;
wght_s = [0.9 0.1];
rate_s = [0.5 4];
[count, label_s] = pm_gen(wght_s, rate_s, T);

% Estimate the parameters using EM
Nit = 10;
[wght,rate,logl,postprob] = pm_em(count, [0.5 0.5], [1 5], Nit);

% Plots
if (S_VER == 1)
  clf;
else
  subplot(2,2,4); % This way you can run the demo more than once in octave!
  clg;
end
% Data
subplot(2,2,1);
if (S_VER == 0) clg; axis; end
title('Data');
plot(count, 'b');
% Log likelihood
subplot(2,2,2);
if (S_VER == 0) clg; axis; end
title('Log-likelihood');
xlabel('iteration');
plot((0:Nit-1), logl, 'b');
% Plot estimated a posteriori probabilities
subplot(2,2,3);
if (S_VER == 0) clg; end
axis([1 T 0 1.1]);
grid;
title('Posterior probabilities');
hold on;
plot(postprob(:,2), 'b');
% Display simulated hidden states
plot((label_s == 2), 'or');

fprintf(1, 'Hit any key to continue.\n');
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '\nPoisson HMM.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
T = 80;
TRANS_s = [0.9 0.1; 0.3 0.7];
rate_s = [0.5 4];
[count, label_s] = ph_gen(TRANS_s, rate_s, T);

% Estimate the parameters using EM
Nit = 10;
[TRANS,rate,logl,postprob] = ph_em(count, ones(2)/2, [1 5], Nit);
% Compute the optimal hidden state sequence (for classification) using Viterbi
class = ph_vit(count, TRANS, rate);

% Plot
if (S_VER == 1) clf; end
% Data plot
subplot(2,2,1);
if (S_VER == 0) clg; axis; end
title('Data');
plot(count, 'b');
% Log-likelihood
subplot(2,2,2);
if (S_VER == 0) clg; axis; end
title('Log-likelihood');
xlabel('iteration');
plot((0:Nit-1), logl, 'b');
% Marginal posterior probabilities
subplot(2,2,3);
if (S_VER == 0) clg; end
axis([1 T 0 1.1]);
grid;
title('Posterior probabilities, MAP sequence');
hold on;
plot(postprob(:,2), 'b');
% Estimated classification
plot(class-1, 'm');
% Display simulated hidden states
plot((label_s == 2), 'or');

fprintf(1, 'Hit any key to continue.\n');
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '\nNegative binomial HMM.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 120;
TRANS_s = [0.9 0.1; 0.3 0.7];
alpha_s = [2 4];
beta_s  = [1 0.25];
[count, label_s] = nbh_gen(TRANS_s, alpha_s, beta_s, T);

% And estimate the parameters using EM
Nit = 10;
[TRANS, alpha, beta, logl, postprob] = nbh_em(count, ones(2)/2, ...
  [1 20], [1 1], Nit);
% Compute the optimal hidden state sequence (for classification) using Viterbi
class = nbh_vit(count, TRANS, alpha, beta);

% Plots
if (S_VER == 1) clf; end
% Data
subplot(2,2,1);
if (S_VER == 0) clg;  axis; end
title('Data');
plot(count,'b');
% Log likelihood
subplot(2,2,2);
if (S_VER == 0) clg;  axis; end
title('Log-likelihood');
xlabel('iteration');
plot((0:Nit-1), logl, 'b');
% Marginal posterior probabilities
subplot(2,2,3);
if (S_VER == 0) clg; end
axis([1 T 0 1.1]);
grid;
title('Posterior probabilities, MAP sequence');
hold on;
plot(postprob(:,2), 'b');
% Estimated classification
plot(class-1,'m');
% Display simulated hidden states
plot((label_s == 2), 'or');

% Plot the marginal distribution (in the stationnary regime)
% Compute negative binomial distributions for all model states
t = 0:max(count);
[tmp1, tmp2, tmp3, tmp4, tmp5, dens] = nbh_em(t, TRANS, alpha, beta, 1);
% and stationary distribution
w = statdis(TRANS);
% Plot estimate of marginal probabilities
subplot(2,2,4);
if (S_VER == 0) clg; axis; end
title('Estimated marginal distribution');
hold on;
plot(t, sum(dens' .* (w(:)*ones(1,length(t)))), '*b');
% Plot empirical estimated probabilities
dhist = zeros(1,length(t));
for i = t
  dhist(1+i) = sum(count == i)/T;
end
plot(t, dhist, '+r');
