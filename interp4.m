function Y = interp4(X, T)
% Y = interp4(X, T)
%   4-point (cubic) interpolation.  X is a sequence; T is a (vector
%   of) real-valued indices on a one-basis.  Y returns the values
%   of X interpolated at each T value.
%   This is based on the lagrange polynomial described on p. 46 of 
%   The Theory and Technique of Electronic Music by Miller Puckette.
% 2014-04-29 Dan Ellis dpwe@ee.columbia.edu

if size(X,1) > 1; X = X'; end

lX = length(X);
Tr = min(lX-3, max(1, floor(T)));

Td = max(0, min(1, T - Tr));

Y = -Td.*(Td-1).*(Td-2)/6.*X(Tr) ...
    + (Td+1).*(Td-1).*(Td-2)/2.*X(Tr+1) ...
    - (Td+1).*Td.*(Td-2)/2.*X(Tr+2) ...
    + (Td+1).*Td.*(Td-1)/6.*X(Tr+3);
