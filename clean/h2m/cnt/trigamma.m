function p =  trigamma(x)

%trigamma  Computes the trigamma function (d^2 log gamma).
%          Use:  p = trigamma(x) where x is a scalar or vector.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 16/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Adapted from GAUSS code by Paul L. Fackler available at
%  http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/loggamma.src

% Formulas 6.4.11 and 6.4.12 with recurrence formula 6.4.6 from
% Abromowitz and Stegun, Dover (1965)

if (any(x <= 0))
  error('trigamma requires positive arguments.');
end

x = x+6;
p = 1./(x.^2);
p = (((((0.075757575757576*p-0.033333333333333).*p+0.0238095238095238) ...
  .*p-0.033333333333333).*p+0.166666666666667).*p+1)./x+0.5*p;
p = 1./((x-1).^2)+1./((x-2).^2)+1./((x-3).^2)+1./((x-4).^2)+...
  1./((x-5).^2)+1./((x-6).^2)+p;
