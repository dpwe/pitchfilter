function y = env_bpf(x, z_r, p_r, w, h)
% y = env_hpf(x, z_r, p_r, w, h)
%    Enhance waveform x by converting into STFT domain using
%    w-point windows hopped at h points, then applying a band-pass
%    filter along time (rasta-like) to each frequency band (in dB
%    domain), then reconstructing.   z_r is zero radius, p_r is
%    pole radius.
% 2014-07-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; z_r = 0.8; end
if nargin < 3; p_r = 0.3; end
if nargin < 4; w = 256; end
if nargin < 5; h = w/4; end

% filter - zeros at dc, f_nyq, one real pole
b = conv([1 -z_r], [1 1]);
a = [1 -p_r];

%freqz(b,a);

X = stft(x, w, w, h);

Xangle = angle(X);
Xdb = 20 * log10(abs(X));
% Floor at 40 dB below peak
marg = 40.0;
Xdb = max(Xdb, repmat(max(Xdb,[],2)-marg, 1, size(Xdb, 2)));


Xhpf = filter_by_row(b, a, Xdb);

y = istft(10.^(Xhpf/20) .* exp(j*Xangle), w, w, h);


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = filter_by_row(B, A, X)
% Apply a single filter to every row of X

nr = size(X,1);
for i = 1:nr
  Y(i,:) = filter(B, A, X(i,:));
end
