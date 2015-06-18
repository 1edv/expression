function dens = c_dgaus(X, mu, Sigma)

%c_dgaus   Computes a set of multivariate normal density values
%	   in the case of diagonal covariance matrices (MEX-file).
%	Use : dens = c_dgaus(X, mu, Sigma) where
%	X (T,p)		T observed vectors of dimension p
%	mu (N,p)	N mean vectors
%	Sigma (N,p)	diagonals of the covariance matrices
%
%	dens (T,N)	density values for each observation and each density.
%
%       c_dgaus is a MATLAB V5 mex-file (to be compiled for your system with
%       the mex command)

% H2M Toolbox, Version 2.0
% Olivier Cappé, 12/03/96 - 06/06/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris
