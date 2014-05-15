function Y = sgram_enhance(X, W_T, W_F)
% Y = sgramenhance(X, W_T, W_F)
%   X is a spectrogram.  Enhance in both time and frequency using a
%   time window of W_T cols and a frequency window of W_F rows.
% 2014-05-14 Dan Ellis dpwe@ee.columbia.edu

% Idea: 
%  - average along time to smooth out noise
%  - local high-pass filter along frequency to suppress local noise

n_fft = 256; % 32 ms at 8 kHz
n_hop = n_fft/4;

XS = stft(X, n_fft, n_fft, n_hop);

XM = abs(XS);

% Normalize values to be ~ 1
scale = sqrt(mean(mean(XM.^2)));
XM = XM/scale;

% median filter along time
XMMFT = medfilt1(XM, W_T, [], 2);

% Enhancement across frequency (in magnitude domain)

if W_F > 1
  % Form local average in freq: 
  % convolution kernel is [1/(n-1) 1/(n-1) .. 0 .. 1/(n-1) 1/(n-1)]
  fkern = ones(1, W_F)/(W_F-1);
  fkern(floor(W_F/2)+1) = 0;
  % Convolve columns to get average of surrounding points
  tkern = [1];
  XFAV = conv2(fkern, tkern, XMMFT, 'same');

  figure(2)
  subplot(211)
  imgsc(XMMFT)
  subplot(212)
  imgsc(XFAV)
  linkaxes();
  figure(1)
  
  
  % Subtract it away (soft max, transition ~ 0.1 units)
  %trans = 0.1;
  %XMLP = trans * log(1+exp((XMMFT - 1.5*XFAV)/trans));
  % hard max
  XMLP = max(0, XMMFT - XFAV);
else
  XMLP = XMMFT;
end
  
% Modify magnitudes
Y = istft(scale*XS.*XMLP./XM, n_fft, n_fft, n_hop);
