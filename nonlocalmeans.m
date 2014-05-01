function DO = nonlocalmeans(D,SR,N,NFFT)
% DO = nonlocalmeans(D,SR,N,NFFT)
%   Rewrite STFT magnitudes using nonlocal means average of STFT
%   frames.  N is how many points to average; NFFT is FFT size.
% 2014-04-28 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; N = 8; end
if nargin < 4; NFFT = 512; end

% Take initial STFT
nhop = NFFT/4;
SD = stft(D, NFFT, NFFT, nhop);
MD = abs(SD);

% Build similarity matrix
% project into log-F domain, distance on cube-root compressed
nlogbin = 72;
fmin = 50;
mx = fft2logfmx(NFFT,SR,nlogbin,fmin);
LD = mx*MD;
expnt = 0.3;
% Simmx stacks 4 successive frames
SM = simmx(abs([LD;rot(LD,1);rot(LD,2);rot(LD,3);rot(LD,4);...
                rot(LD,5);rot(LD,6);rot(LD,7)]).^expnt);

imgsc(SM)

% Find nearest neighbors of each frame
[vv,xx] = sort(SM, 'descend');

% Build up average of nearest neighbors
AD = zeros(size(SD,1), size(SD,2));
for i = 1:N; 
  AD = AD + MD(:, xx(i+1,:)); 
end
AD = AD/N;

% Reconstruct using original phase
DO = istft(AD.*angle(SD), NFFT, NFFT, nhop);
