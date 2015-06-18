function w = hanning(n)

%hanning   Hanning window (introduced for MATLAB compatibility).
%          Use: w = hanning(n).
%          Note: As of octave 2.0.14, this function is not needed
%          anymore (there is an octave function of the same name).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 15/04/99 - 17/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

w = sin(pi*(0:n-1)'/n).^2;
