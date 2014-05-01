function y = pitchfilter(d, sr)
% Y = pitchfilter(D,SR)
%   Enhance a signal by filtering at the harmonics of the detected
%   pitch.  D@SR is a waveform; run noise-robust pitch tracking
%   (via SAcC), then resample the waveform so each frame maps the
%   detected pitch to 100 Hz.  Then apply a comb filter to
%   emphasize all the harmonics of 100 Hz, then undo the resampling
%   to get back to the original pitch.
% 2014-05-01 Dan Ellis dpwe@ee.columbia.edu

doplot = 1;

% (1) Run SAcC pitch tracker with the noisy rats classifier, and
% with the unvoiced state discounted, to get a near-continuous
% pitch track.
SAcCdir = '/u/drspeech/data/RATS/code/SAcC';
if exist('SAcC') == 0
  addpath(SAcCdir);
end
configfile = 'conf/rats_sr8k_bpo6_sb24_k10.config';
params = config_read_srs(fullfile(SAcCdir, configfile));
% fix up paths
for param = {'pca_file', 'wgt_file', 'norms_file', 'pcf_file'}
  val = getfield(params, param{1});
  params = setfield(params, param{1}, fullfile(SAcCdir, val));
end
% modify voicing prior
params.hmm_vp = 0.01;

% Run SAcC pitch tracker
[pitch, lobs, pvx, times] = SAcC_main(d, sr, params);
delta_t = times(2) - times(1); % Assume it's fixed throughout
% Better alignment of pitch to signal?
times = times + 0.015;

% Median filter the prob(voicing) for smoother crossfades
pvxwin = 0.2;
% Force window to be an odd number of points
pvxwinpts = round(pvxwin/(2*delta_t)) * 2 + 1;
pvx = medfilt1(pvx, pvxwinpts);

if doplot
  figure(1)
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
  axis([0 length(d)/sr 0 2000]);
end

% Build the time mapping
% Map is a two rows; first row indicates a point in original space
% corresponding element in second row is where it ends up
tt0 = times(1);
target_pitch = 100.0;
vmap = [times; zeros(1, length(times))];
for i = 1:length(pitch)
  vmap(2,i) = tt0; 
  % map any unvoiced frames to neutral resampling (100 Hz)
  tt0 = tt0 + (pitch(i) + target_pitch*(pitch(i)==0))/target_pitch ...
               * delta_t;
end

% resample waveform to 100 Hz
tt = [1:length(d)]/sr;
% figure the value of the original waveform at the inverse-mapped
% times corresponding to regular samples in the mapped waveform
% times of samples in modified domain
ttd = 0:(1/sr):max(vmap(2,:));
% map these back to sample times in the original
ttm = map_vals(ttd, inv_map(vmap));
% interpolate the base waveform to get those values
dm = interp4(d', ttm*sr - 2);  % -2 fixes for the 1 sample delay
                               % introduced in each interp4

if doplot
  subplot(412)
  specgram(dm, sg_fft, sr, sg_fft, sg_olp);
  title('resampled to pitch = 100 Hz');
  caxis(cax);
  colormap(1-gray);
  axis([0 length(dm)/sr 0 2000]);
end

%y = dm;

% Apply 100 Hz enhancement
% Simple comb - feedback at the period, maybe a zero nearby to
% limit impact
f_ord = round(sr/target_pitch);
f_a = zeros(1, f_ord+1);
pole_radius = 0.995;
zero_radius = 0.1;
f_a(1) = 1;
%f_a(f_ord+1) = -(pole_radius^(f_ord));
f_b(1) = 1;
%f_b(f_ord+1) = -(zero_radius^(f_ord));
f_b(f_ord+1) = 1;
% Single FIR
fir_ord = 16;
cutoff = 0.2;
b = fir1(fir_ord, cutoff);
% Comb-ize
f_b = zeros(1, 1 + length(b)*f_ord);
for i = 1:length(b)
  f_b(1 + (i-1)*f_ord) = b(i);
end
f_a = 1;

% Apply filter
%dmf = filter(f_b, f_a, dm);
dmf = conv(f_b, dm);

% Remove linear phase delay
dmf = dmf(((length(b)-1)/2*f_ord)+[1:length(dm)]);

if doplot
  subplot(413)
  specgram(dmf, sg_fft, sr, sg_fft, sg_olp);
  title('Filtered to enhance 100 Hz');
  caxis(cax);
  colormap(1-gray);
  axis([0 length(dm)/sr 0 2000]);
end

% Resample back to original domain
ttu = map_vals(tt, vmap);
du = interp4(dmf, ttu*sr);

% crossfade based on pvx
dsfact = round(sr*delta_t);
pvxi = resample(pvx, dsfact, 1);

maxl = max([length(d), length(du), length(pvxi)]);
d(maxl) = 0;
du(maxl) = 0;
pvxi(maxl) = 0;

du = pvxi.*du' + (1-pvxi).*d;

if doplot
  subplot(414)
  specgram(du, sg_fft, sr, sg_fft, sg_olp);
  title('Resampled back to original pitch and mixed with original');
  caxis(cax);
  colormap(1-gray);
  axis([0 length(du)/sr 0 2000]);

  figure(2)
  freqz(f_b, f_a);
  %zplane(f_b, f_a);
  
end

y = du;
