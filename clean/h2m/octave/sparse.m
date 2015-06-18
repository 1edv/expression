function y = sparse(x)

%sparse    Does nothing (introduced for MATLAB compatibility).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 15/04/99 - 14/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Note: This will badly fail if not using the --traditional switch
global FAKE_SPARSE_WARNING_MESSAGE_ALREADY_PRINTED_ONCE;

if isempty(FAKE_SPARSE_WARNING_MESSAGE_ALREADY_PRINTED_ONCE)
  fprintf(1,'Warning: I don''t know how to handle sparse matrices, matrix will be full\n');
  FAKE_SPARSE_WARNING_MESSAGE_ALREADY_PRINTED_ONCE = 1;
end
y = x;
