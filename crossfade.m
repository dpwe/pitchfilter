function y = crossfade(a, b, g, gdsfact)
% y = crossfade(a, b, g, gdsfact)
%   Y returns as a crossfade between signals a and b according to 
%   mix factor g : y = g*a + (1-g)*b.  But g is at a sampling rate 
%   gdsfact lower than a and b; also take care of making everything
%   the same length (max(length(a), length(b))).
% 2014-05-14 Dan Ellis dpwe@ee.columbia

%dsfact = round(sr*median(diff(times)));

gi = resample(g, gdsfact, 1);

maxl = max([length(a), length(b)]);

a(maxl+1) = 0;
b(maxl+1) = 0;
gi(maxl+1) = 0;
gi = gi(1:maxl+1);

y = gi.*a + (1-gi).*b;

y = y(1:maxl);
