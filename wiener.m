function Y = wiener(X, SR)
% Y = wiener(X, SR)
%      Apply wiener filter enhancement, attempting to duplicate
%      ICSI's "nr" process.  
% 2014-05-15 Dan Ellis dpwe@ee.columbia.edu

% preemphasize
emphfilt = [1 -0.96];
X = filter(emphfilt, 1, X);

% STFT
targetwinsec = 0.025;
nfft = 2^round(log(targetwinsec*SR)/log(2));
nhop = nfft/4;
fftframesec = nhop/SR;

% Initial spectrum
XS = stft(X, nfft, nfft, nhop);

% Magnitude
XMag = abs(XS);

% Mel domain
nmel = 40;
WW = fft2melmx(nfft, SR, nmel, 1, 0, SR/2, 0, 1);
XMel = WW * XMag;

% Figure voice activity (simple smoothed energy threshold)
VAD = yet_another_vad(XMel, SR/nhop);

% "Modified Wiener filter" per icslp02-aurora.pdf 
% "Qualcomm-ICSI-OGI features for ASR", Adami et al, Interspeech 2002.

% Noise spectrum estimate

% Guess noise from first <initialestwin> frames
%initialestwin = 100;
%noiseframes = 1:initialestwin;
%
% No, estimate from all non-VAD frames
noiseframes = find(VAD);

% Estimate average noise over all noiseframes
% In log domain, in keeping with auroralib.c
Noise = idB(mean(dB(XMel(:, noiseframes)),2));

% Duplicate for all time to give noise estimate What
What = repmat(Noise, 1, size(XMel,2));

% Online estimate of noise floor
%pole_r = 0.98;
%What = filter_by_row((1 - pole_r), [1 -pole_r], XMel, pole_r*Noise);

X2 = XMel.^2;
What2 = What.^2;

% Wiener filter

% SNRapost is actually (signal+noise)/(noise estimate) (noiscomp.c:138)
%SNRapost = 10*log10(max(1e-2, (sum(X2) - sum(What2))./sum(What2)));
SNRapost = 10*log10(max(1e-2, (sum(X2))./sum(What2)));

% eqn (2) from paper
gamma_k = max(1.25, min(3.125, -1.875/20*SNRapost + 3.125));

% eqn (1) from paper
beta = 0.01;
Hinst2 = max(beta, (X2 - repmat(gamma_k, nmel, 1) .* What2)./X2);

% Smooth in time and frequency
twinlen = 21;
t_kern = hanning(twinlen)'/sum(hanning(twinlen));
fwinlen = 5;
f_kern = hanning(fwinlen)'/sum(hanning(fwinlen));

% nr.c actually uses a one-pole smoothing along time with 
% alpha = 0.1, and a two-frame advance to account for lag. 
% I think it uses 16 ms hop.
% (noisecomp.c:164; nr.c:251; 
% freqfilt looks to be boxcar smoothing over 21 FFT bins
% (noisecomp.c:238 et seq.)

H2 = conv2(f_kern, t_kern, Hinst2, 'same');

alpha = 0.001;

Shat = sqrt(max(alpha*What2, X2 .* H2));

% Back into FFT domain


% Reconstruct
XMask = WW'*(Shat./XMel);

Y = istft(XS.*XMask, nfft, nfft, nhop);

% deemphasize
Y = filter(1, emphfilt, Y);


% plot
%specgram(X, nfft, SR);
ax = [-30 30];

subplot(411)
imgsc(dB(XMel));
caxis(ax);
title('X')

subplot(412)
imgsc(dB(What));
caxis(ax);
title('What')

subplot(413)
%plot(SNRapost);
%title('SNRapost');
imgsc(dB(Hinst2)/2)
caxis(ax);
title('Hinst');

subplot(414)
%imgsc(dB(Shat));
%caxis(ax);
%title('Shat');
imgsc(dB(H2)/2)
caxis(ax);
title('H');

linkaxes

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = filter_by_row(B, A, X, Z)
% Apply a single filter to every row of X; use initial state

nr = size(X,1);
for i = 1:nr
  Y(i,:) = filter(B, A, X(i,:), Z(i,:));
end
