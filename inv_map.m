function Y = inv_map(M)
% function Y = inv_map(M)
%    construct the inverse of a monotonic map s.t. 
%    map_vals(map_vals(vals, vmap), inv_map(vmap)) == vals
% 2014-04-29 Dan Ellis dpwe@ee.columbia.edu

% you simply interchange the input points and the output points
Y = M([2 1], :);
