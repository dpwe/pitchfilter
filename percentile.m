function v = percentile(d,n)
% v = percentile(d,n)
%    Return for each col of v the n'th percentile where 0<n<1
% 2004-10-04 dpwe@ee.columbia.edu

[nr,nc] = size(d);

if nr == 1
  d = d';
  nr = nc;
  nc = 1;
end

x = sort(d);
v = x(1+floor(n*nr),:);
