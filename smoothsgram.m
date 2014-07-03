function Y = smoothsgram(X, SR)
% Y = smoothsgram(X, SR)
%    Reduce noise in STFT domain by smoothing along time
% 2014-05-22 Dan Ellis dpwe@ee.columbia.edu

nfft = 256;
% Use rectangular window because we're expecting bin-aligned harmonics
%win = ones(1,nfft)/nfft;
%nhop = nfft/2;
win = 256;
nhop = nfft/4;

D = stft(X, nfft, win, nhop);
[nrows, ncols] = size(D);

old = 0;
if old
  % old approach based on dividing & predicting (Bello-style)
  
  % Calculate the actual phase advance for successive windows
  BonA = D(:, 2:end)./D(:, 1:end-1);

  % Prediction of next frame by extrapolating previous two
  Dpred = zeros(nrows, ncols);
  Dpred(:, 3:end) = D(:, 2:end-1) .* BonA(:, 1:end-1);
  % Plot comparisons
  cax = [-30 30];

  subplot(411)
  imgsc(dB(D))
  caxis(cax)
  title('D original');

  subplot(412)
  imgsc(dB(Dpred))
  caxis(cax)
  title('D predicted');

  subplot(413)
  imgsc(dB(D - Dpred))
  caxis(cax)
  title('Error D - Dpred');

  subplot(414)
  imgsc(dB(D - Dpred) - dB(D))
  caxis(cax)
  title(' (D - Dpred)/D ')

end

% Half-window for smoothing
P = 5;
% Local weighting window
W = hann(2*P+1)';
W = W ./ sum(W);

% Calculate phase advance by multiplying adjacent conjugates
% then small magnitudes are downweighted
phaseadv = D .* conj(D(:, [1, 1:end-1]));
% Discount some of the magnitude weighting
phaseadv = phaseadv ./ sqrt(abs(phaseadv));

% Now smooth
smphaseadv = conv2(W, [1], phaseadv, 'same');
% Angle of this is the smoothed local phase advance
alpha = smphaseadv./abs(smphaseadv);

% Each complex STFT bin is the average of a neighborhood in 
% time, rotated to the current bin by a smoothed average phase
% increment

Dpred = zeros(nrows, ncols);

for j = (P+1):(ncols-P)
  Dpred(:, j) = sum(bsxfun(@times, W, D(:, [(j-P):(j+P)])) ... 
                    ./ bsxfun(@power, alpha(:, j), [-P:P]), 2);
end

% Reconstruct with hann windows
Y = istft(Dpred, nfft, win, nhop);