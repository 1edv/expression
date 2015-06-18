function [count, label] = ph_gen(TRANS, rate, T)

%ph_gen	  Simulates data from a Poisson HMM.
%         Use: [count,label] = ph_gen(TRANS,rate,T|label) where count
%         is an array of length T which contains the simulated data
%         and label contains the corresponding simulated state
%         sequence. If the last argument is of length one, it is taken
%         as the number of observations, otherwise it is considered as
%         a specified state sequence.
%         Requires Matlab's Statistics toolbox or GNU Octave.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 30/12/97 - 22/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (~((exist('randindx') == 2) | (exist('randindx') == 3)))
  error('Function randindx (from h2m toolbox) is missing.');
end
% Check that we have something to simulate from Poisson
if (exist('poisson_rnd') == 2)
  S_VER = 0; % Octave
elseif (exist('poissrnd') == 2)
  S_VER = 1; % Matlab + statistics toolbox
else
  error(['Either poissrnd (Matlab''s Statistics toolbox) or\n' ...
         'poisson_rnd (Octave >= 2.0.14) are required (but read\n' ...
	 'the H2M doc to see how to get around the problem).']);
end

% Inputs arguments
error(nargchk(3, 3, nargin));
N = ph_chk(TRANS, rate);
% Make sure that the output will be a column vector
rate = reshape(rate, N, 1);

if (length(T) == 1)
  % T contains the length of data to simulate
  % First simulate labels
  label = zeros(T,1);
  % Simulate initial state
  label(1) = randindx(ones(1,N)/N, 1, 1);
  % Use Markov property for the following time index
  for t = 2:T
    label(t) = randindx(TRANS(label(t-1), :), 1, 1);
  end
else
  % T directly contains a sequence of labels
  label = reshape(T, length(T), 1);
end

% Draw the Poisson data
if (S_VER == 0)
  count = poisson_rnd(rate(label));
else
  count = poissrnd(rate(label));
end
