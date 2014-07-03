function y = env_hpf(x, r, w, h)
% y = env_hpf(x, r, w, h)
%    Enhance waveform x by converting into STFT domain using
%    w-point windows hopped at h points, then applying a high-pass
%    filter along time [1 -r] to each frequency band (in dB
%    domain), then reconstructing.
% 2014-07-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; r = 0.5; end
if nargin < 3; w = 256; end
if nargin < 4; h = w/4; end

X = stft(x, w, w, h);

Xdb = 20 * log10(abs(X));
Xangle = angle(X);

Xhpf = filter_by_row([1 -r], [1], Xdb);

y = istft(10.^(Xhpf/20) .* exp(j*Xangle), w, w, h);


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = filter_by_row(B, A, X)
% Apply a single filter to every row of X

nr = size(X,1);
for i = 1:nr
  Y(i,:) = filter(B, A, X(i,:));
end
