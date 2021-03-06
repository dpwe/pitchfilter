function Y = wiener_icsi(X, SR)
% Y = wiener_icsi(X, SR)
%      Apply wiener filter enhancement, attempting to duplicate
%      ICSI's "nr" process.  
%      This version is as close as possible to ICSI's
% 2014-05-15 Dan Ellis dpwe@ee.columbia.edu

% STFT
targetwinsec = 0.025;
nfft = 2^round(log(targetwinsec*SR)/log(2));
nhop = nfft/2;
fftframesec = nhop/SR;

XS = stft(X, nfft, nfft, nhop);
% Magnitude
XMag = abs(XS);

% No Mel in ICSI processing

% Figure voice activity (simple smoothed energy threshold)
VAD = yet_another_vad(XMag, SR/nhop);

% "Modified Wiener filter" per icslp02-aurora.pdf 
% "Qualcomm-ICSI-OGI features for ASR", Adami et al, Interspeech 2002.

% Noise spectrum estimate - simple (dB) average over all non-VAD frames
% (log domain to match auroralib.c)
Noise = exp(mean(log(XMag(:, find(VAD))),2));

% Duplicate for all time to give noise estimate W_hat
W_hat = repmat(Noise, 1, size(XMag, 2));

X2 = XMag.^2;
W_hat2 = W_hat.^2;

% Wiener filter

% SNRapost is actually (signal+noise)/(noise estimate) (noiscomp.c:138)
SNRapost = 10*log10(max(1e-2, (sum(X2))./sum(W_hat2)));

% Mapping SNR to overmasking factor (eqn (2) from paper)
% (over_est_factor at noiscomp.c:149)
gamma_k = max(1.25, min(3.125, -1.875/20*SNRapost + 3.125));

% Estimating masking filter as overmasked SNR (eqn (1) from paper)
% (noisecomp.c:160)
beta = 0.01;
Hinst2 = max(beta, (X2 - repmat(gamma_k, size(X2, 1), 1) .* W_hat2)./X2);
% This gets squared at noisecomp.c:162!
%Hinst2 = Hinst2 .^ 2;
% Sounds better without it - dips too severe if present

% Smooth in time (noisecomp.c:164)
% Original values (time_delay needs to be longer for smaller alpha)
%t_alpha = 0.1;
%time_delay = 2;
% These values sound better to me (for BABEL_OP1_204_MIC)
t_alpha = 0.4;
time_delay = 1;
% Simple 1-pole filter on each row
H_smoo = filter_by_row(t_alpha, [1 -(1-t_alpha)], Hinst2);

% Smooth in frequency
fwinlen = 21;  % 21 in original
% Original is boxcar (noisecomp.c:238 et seq.)
%f_kern = ones(fwinlen, 1)/fwinlen;
% Tapered window is effectively narrower, sounds better
f_kern = hann(fwinlen)/sum(hann(fwinlen));
H2 = conv2(f_kern, [1], H_smoo, 'same');

% Time-advance mask (e.g. noisecomp.c:94)
H2 = [H2(:, (time_delay+1):end), repmat(H2(:, end), 1, time_delay)];

% Actual filtering
% worst-case minimum (noiscomp.c:288, nr.c:255)
alpha = 0.001;
Shat = sqrt(max(alpha*W_hat2, X2.*H2));

% Back out of STFT domain
Y = istft(XS.*Shat./XMag, nfft, 0, nhop);

% Y is my final answer

% Plotting all at end to preserve flow
do_plot = 0;
if do_plot
  % plot
  %specgram(X, nfft, SR);
  ax = [-30 30];

  subplot(411)
  imgsc(dB(XS));
  caxis(ax);
  title('X')

  subplot(412)
  imgsc(dB(Hinst2)/2)
  caxis(ax);
  title('Hinst');

  subplot(413)
  imgsc(dB(H2)/2)
  caxis(ax);
  title('H');

  subplot(414)
  imgsc(dB(Shat));
  caxis(ax);
  title('Shat');

  linkaxes

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = filter_by_row(B, A, X)
% Apply a single filter to every row of X; use initial state

nr = size(X,1);
for i = 1:nr
  Y(i,:) = filter(B, A, X(i,:));
end
