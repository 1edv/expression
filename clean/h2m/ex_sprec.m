%ex_sprec  (Unrealistic) example of isolated word recognition
%          (in speaker-dependent mode) with the h2m functions.
%
%	Beware : This is not a real demo, so you should take a look at the
%       source file before executing it if you wan't to understand what's
%       going on!

% H2M Toolbox, Version 2.0
% Olivier Cappé, 09/04/98 - 17/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

%clear;

% Try to figure out which version of matlab or octave is used (this is crude)
s_ver = version;
if (s_ver(1) == '2')
  S_VER = 0; % Octave 2.0 something
  fprintf(1,['I am assuming that your are using GNU Octave version\n' ...
    '2.0.something in --traditional mode (congratulations!).\n']);
elseif (s_ver(1) == '4')
  S_VER = 1; % Matlab V 4
  fprintf(1,'I am assuming that your are using Matlab, version 4.\n');
else
  S_VER = 2;
  fprintf(1,['I am assuming that your are using Matlab, version 5' ...
  ' (or above).\n']);  
end
tmp = input(['If you are not happy with that enter one of\n  [0] GNU' ...
  ' Octave --traditional\n  [1] Matlab 4\n  [2] Matlab 5 or above\nor just' ...
  ' hit <return> if you believe that the current choice is right : ']);
if (~isempty(tmp))
  S_VER = tmp;
end
if (S_VER == 0)
  path(path, 'octave');
end

% Data base
NREP = 12;      % Number of utterances per word
% Limits of each word in the signal file (in samples)
SEG_FILE = 'data/digits.mat';
% Signal file (in little-endian short format)
SIG_FILE = 'data/digits.sig';

% Analysis parameters
lfen = 256;             % Window size
dec = 128;              % Hop-size
taper = hanning(lfen)'; % Weighting window
pnorm = sum(taper.^2);
p = 10;                 % Cepstrum order

% Training parameters
Ntrain = 4;     % Number of training utterances (the remaining NREP-Ntrain
                % occurences will be used for testing)
DIAG_COV = 1;   % Force use diagonal covariance matrices
QUIET = 1;      % Make training routines silent
N = 5;          % Number of states per word model
% Transition matrix is more or less arbitrary and will not be estimated (too
% few utterances are available)
A = sparse(0.85*diag(ones(1,N))+0.15*diag(ones(1,N-1),1));
A(N,N) = 1;
NIT = 10;       % Number of EM iterations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data file management (take care of problems occuring with big-endian sys.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'Checking data file...\n');
if (S_VER ~= 1)
  % You are using MATLAB 5 or above of GNU octave, no problem...
  % open the signal file
  fid = fopen(SIG_FILE, 'r', 'ieee-le');
else
  % Try to read the first sample
  fid = fopen(SIG_FILE, 'r');
  ts = fread(fid, 1, 'short');
  fclose(fid);
  if (ts ~= -4) % Magic number = value of first sample
    input(['File format is not readable. I can fix the problem by\n' ...
           'overwritting the data file. Abort with ^C now if you\n'...
           'are not happy with this or type <return> to proceed']);
    fid = fopen(SIG_FILE, 'r');
    s = fread(fid, inf, 'char');
    fclose(fid);
    s = reshape([s(2:2:length(s))'; s(1:2:length(s)-1)'], length(s), 1);
    fid = fopen(SIG_FILE, 'w');
    fwrite(fid, s, 'char');
    fclose(fid);
    fprintf(1, 'Data file has been modified\n');
  end
  % Everything should now be correct in both cases, open the file
  fid = fopen(SIG_FILE, 'r');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterization (with standard power cepstrum coeffs.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'Computing short-time parameters...\n');
% Number of samples in data file
status = fseek(fid, 0, 'eof');
l = ftell(fid)/2;
% Number of short-time frames
nbf = 1+fix((l-lfen)/dec);
% Array for cepstrum coefficients
C = zeros(nbf, p);
% Compute standard cepstral coeffs.
x = zeros(1, lfen);
% Read first portion of frame
status = fseek(fid, 0, 'bof');
x(1+dec:lfen) = fread(fid, lfen-dec, 'short');
for n = 1:nbf
  % Read input signal
  x(1:lfen-dec) = x(1+dec:lfen);
  x(lfen-dec+1:lfen) = fread(fid, dec, 'short');
  % Periodogram
  S = abs(fft(x.*taper)).^2/pnorm;
  % Power cepstrum (use normalized transform here)
  S = sqrt(lfen)*ifft(log(S+realmin));
  C(n,:) = real(S(1:p));        % Including c(0)
end
fclose(fid);    % In real life, you would probably save the cepstrum
                % parameters to a file about here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load segmentation file (this is a mat file with indexes of the first an last
% samples corresponding to each word)
load(SEG_FILE);
% For each word
fprintf(1, 'Training model for word:');
for w = 0:9
  fprintf(1, ' %d', w);
  % Training of words w
  % Gather training parameters for the word in a single matrix
  X = [];
  st = [];
  for i = (1+w+((1:Ntrain)-1)*10)
    if (segment(i,1) ~= w)
      error('Something''s wrong with the segmentation file');
    end
    st = [st; size(X,1)+1];
    n = 1+round(segment(i,2:3)/dec);
    X = [X; C(n(1):n(2),:)];
  end
  T = size(X,1);

  % Training. Because there is very few available training data, the covariance
  % matrices are diagonal and shared by all states of the word model (ie. all
  % states have the same covariance matrix given by Sigma)
  [mu,Sigma] = hmm_mint(X, st, N, DIAG_COV,QUIET);
  Sigma = ones(N,1)*mean(Sigma);                        % Shares covariance
  logl = zeros(1, NIT);
  for n = 1:NIT
    % Expectation step of the EM algorithm
    [tmp, logl(n), gamma] = hmm_mest(X, st, A, mu, Sigma, QUIET);
    % Unconstrained Maximization step of the EM algorithm
    [mu, Sigma] = mix_par(X, gamma, DIAG_COV, QUIET);
    % Modification of the EM Maximization due to the constraint that all
    % covariance matrices are identical
    Sigma = ones(N,1)*(sum((sum(gamma)'*ones(1,p)).*Sigma)/T);
  end
  sigma = Sigma(1,:);
  % Write a separate mat file for each words model in subdirectory data
  eval(['save data/mod' int2str(w) ' A mu sigma']);
end
fprintf(1,'\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of tests
Ntest = NREP - Ntrain;
% Matrix of recognized words
reco = zeros(Ntest,10);
fprintf(1, 'Test utter.:\n');
% For each test...
for i = 1:Ntest
  fprintf(1, ' %d', i);
  % ...test all available word
  fprintf(1, ' (word:');
  for t = 0:9
     fprintf(1, ' %d', t);
    if (segment(1+(Ntrain*10)+(i-1)*10+t,1) ~= t)
      error('Something''s wrong with the segmentation file (test)');
    end
    n = 1+round(segment(1+(Ntrain*10)+(i-1)*10+t,2:3)/dec);
    score = zeros(1,10);
    for w=0:9
     if (S_VER == 0)
       eval(['load -force data/mod' int2str(w)]); % Octave call
     else
       eval(['load data/mod' int2str(w)]);
     end
     Sigma = ones(N,1)*sigma;
     score(1+w) = hmm_vit(C(n(1):n(2),:), A, [1 zeros(1,N-1)], mu, Sigma, ...
       QUIET);
     end
  [tmp, reco(i,1+t)] = max(score);
  end
  fprintf(1,')\n');
end
% Print recognized words for all tests (results are in matrix reco)
fprintf(1, 'Recognized words:\n');
for i = 1:length(reco(:,1))
  fprintf(1, 'Utter. %d:', i);
  fprintf(1, ' %d', reco(i,:)-1);
  d = reco(i,:)-(1:10);
  if (any(d))
    fprintf(1, ' - %d error(s)\n', length(find(d ~= 0)));
  else
    fprintf(1, ' - correct\n');
  end
end

fprintf(1,['\nWord models are in directory data: mod0 to mod9.mat.\n' ...
	   'You may erase them now (or they will be overwritten\n' ...
	   ' next time you run this script).\n']);
