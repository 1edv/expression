%digits    Script for plotting and listening to the data in file digits.sig.

% H2M Toolbox, Version 2.0
% Olivier Cappé, 07/02/97 - 15/04/99
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Constants pertaining to C 'short' files
SIZE_OF_SHORT = 2;
MAX_SHORT = 32768;

% Load segmentation file
load digits;

% Load the portion of signal corresponding to each word and listen to it
matlab_version = version;
if (matlab_version(1) == '5')
  % You are using MATLAB 5, no problem...
  fid = fopen('digits.sig', 'r', 'ieee-le');
else
  % Open anyway (will be unpleasant to hear if the machine is not
  % little-endian) !
  fid = fopen('digits.sig', 'r');
end
for i = 1:length(segment(:,1))
  fseek(fid, SIZE_OF_SHORT*segment(i,2), 'bof');
  s = fread(fid, segment(i,3)-segment(i,2)+1, 'short');
  plot(s);
  title(int2str(segment(i,1)));
  drawnow;
  sound(s/MAX_SHORT);
  %pause
end;
fclose(fid);
