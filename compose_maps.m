function Y = compose_maps(M1, M2)
% Y = compose_maps(M1, M2)
%    return a single [x,y] map that represent y = map1(map2(x))
% 2014-04-29 Dan Ellis dpwe@ee.columbia.edu

% Break points will be all the edges in map2, and all the edges in map1 
% when projected through the inverse of map2
mapped_edges = map_vals(M1(1,:), inv_map(M2));
all_edges = unique([M2(1,:), mapped_edges]);
all_vals = map_vals(map_vals(all_edges, M2), M1)
Y = [all_edges; all_vals];
