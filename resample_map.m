function y = resample_map(x, sr, vmap)
% y = resample_map(x, sr, vmap)
%    Interpolate x to find y s.t. times in y relate to times in x
%    according to the piecewise-linear map defined in vmap 
%    (times in seconds; first row is times in x, second row is
%    corresponding times in y).
% 2014-05-14 Dan Ellis dpwe@ee.columbia.edu

xdur = length(x)/sr;
% figure the value of the original waveform at the inverse-mapped
% times of samples in modified domain, which are...
ttd = 0:(1/sr):map_vals(xdur, vmap);
% map these back to sample times in the original
ttm = map_vals(ttd, inv_map(vmap));
% interpolate the base waveform to get those values
y = interp4(x, ttm*sr);
