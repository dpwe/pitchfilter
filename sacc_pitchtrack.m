function [pitch, pvx, times] = sacc_pitchtrack(d,sr)
% [pch, pvx, times] = sacc_pitchtrack(d,sr)
%    Run the SAcC pitch tracker on noisy utterance.
%    pch returns pitch in Hz, pvx returns prob(voicing), 
%    times is sample times.
% 2014-05-14 Dan Ellis dpwe@ee.columbia.edu

doplot = 1;

% Add SAcC functions to path
SAcCdir = '/u/drspeech/data/RATS/code/SAcC';
if exist('SAcC') == 0
  addpath(SAcCdir);
end

% Load the high-noise configuration parameters
configfile = 'conf/rats_sr8k_bpo6_sb24_k10.config';
params = config_read_srs(fullfile(SAcCdir, configfile));

% fix up paths
for param = {'pca_file', 'wgt_file', 'norms_file', 'pcf_file'}
  val = getfield(params, param{1});
  params = setfield(params, param{1}, fullfile(SAcCdir, val));
end

% modify voicing prior to minimize number of unvoiced frames
params.hmm_vp = 0.001;  %0.01;

% Run SAcC pitch tracker
[pitch, lobs, pvx, times] = SAcC_main(d, sr, params);
times = times';
dt = times(2) - times(1); % Assume it's fixed throughout
% Better alignment of pitch to signal?
times = times + 0.015;

% Median filter the prob(voicing) for smoother crossfades
pvxwin = 0.2;  % seconds
% Force window to be an odd number of points
pvxwinpts = round(pvxwin/(2*dt)) * 2 + 1;
pvx = medfilt1(pvx, pvxwinpts);

if doplot
  %figure(1)
  subplot(411)
  sg_fft = 1024;
  sg_olp = sg_fft - sg_fft/8;
  specgram(d, sg_fft, sr, sg_fft, sg_olp);
  title('Noisy signal');
  cax = [-60 20];
  caxis(cax);
  colormap(1-gray);
  hold on;
  plot(times, pitch, '-r', times, pvx*1000, '-g');
  legend('pitch', 'smoothed pvx');
  hold off;
  % just pitch region
  axis([0 length(d)/sr 0 1000]);
end
