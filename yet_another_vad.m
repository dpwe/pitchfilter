function VAD = yet_another_vad(M, FR, SW)
% V = yet_another_vad(M, FR, SW)
%    M is a spectrogram-like (linear magnitude) array, with columns
%    at frame rate FR.  Estimate an energy-based VAD, return in V a
%    binary mask indicating voice activity with one entry per
%    column. 
%    SW is the smoothing window in seconds (default 0.2).
% 2014-05-15 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  SW = 0.2;  end

% Figure energy and use for VAD
E = sqrt(sum(M.^2));
smwinfrm = round(SW*FR);

Esm = conv2([1], hann(smwinfrm)/sum(hann(smwinfrm)), E, 'same');

Ethr = 1.5*percentile(Esm, 0.1);
VAD = Esm > Ethr;

% Plot VAD diagnostics
if 0
  subplot(411);
  imgsc(dB(M)); caxis([-30 30]);
  ax = axis();
  subplot(412);
  plot(1:length(Esm), Esm, 1:length(VAD), 10*VAD, '-r')
  axis([ax(1), ax(2), min(Esm), max(Esm)]);
  subplot(413);
  hist(Esm, 100);
end
