function y = table1(T,x)

%table1    Linear interpolation from table (introduced for MATLAB compatibility).
%          Use: y = table1(T,x).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 31/12/98 - 31/12/98
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

y = zeros(size(x));

for i = 1:length(x)
  xi = x(i);
  if (any(T(:,1) == xi))
    y(i) = T(find(T(:,1) == xi),2);
  else
    tinf = find(T(:,1) < xi);
    tinf = tinf(length(tinf));
    tsup = find(T(:,1) > xi);
    tsup = tsup(1);
    alpha = (xi-T(tinf,1))/(T(tsup,1)-T(tinf,1));
    y(i) = T(tinf,2)*(1-alpha) + T(tsup,2)*alpha;
  end
end
