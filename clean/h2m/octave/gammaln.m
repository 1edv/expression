function z = gammaln(x)

%gammaln   Computes the log of the gamma function.
%          Use:  z = gammaln(x) where x is a scalar or vector.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 17/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Adapted from GAUSS code by Paul L. Fackler available at
%  http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/loggamma.src

% Original comments follow:
% Source: Pike, M.C., and I.D. Hill, Algorithm 291, Communications
% of the ACM, 9,9:p.684 (Sept, 1966). Accepts X>0, accuracy to 10
% decimal places.  Uses Sterling's formula.

if (any(x <= 0))
  error('gammaln requires positive arguments.');
end

x = x+6;
z = 1./(x.^2);
z = (((-0.000595238095238*z+0.000793650793651) ...
  .*z-0.002777777777778).*z+0.083333333333333)./x;
z = (x-0.5).*log(x)-x+0.918938533204673+z ...
  -log(x-1)-log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
