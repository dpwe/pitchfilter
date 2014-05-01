function Y = map_vals(X,vmap)
% Y = map_vals(X,vmap)
%     X is a vector, map is two rows indicating input, output pairs 
% 2014-04-29 Dan Ellis dpwe@ee.columbia.edu

% Just run on a vector
% figure slopes
gap = diff(vmap, 1, 2);
% repeat final gap on both src and dst rows
gap = [gap, gap(:,end)];
% Don't allow zero (or negative?) gaps
gap(find(gap <= 0)) = eps;
% do mapping.  1+ for Matlab indexing later on
Xix = 1 + max(0, sum(bsxfun(@gt, X', vmap(1,:)), 2)-1);
Xdelta = (X - vmap(1, Xix))./gap(1, Xix);
Y = vmap(2, Xix) + Xdelta.*gap(2, Xix);

