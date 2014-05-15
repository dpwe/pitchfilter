function [y,f,g,m] = pitchfilter(d, sr)
% [Y,F,G,T,M] = pitchfilter(D,SR)
%   Enhance a signal by filtering at the harmonics of the detected
%   pitch.  D@SR is a waveform; run noise-robust pitch tracking
%   (via SAcC), then resample the waveform so each frame maps the
%   detected pitch to 100 Hz.  Then apply a comb filter to
%   emphasize all the harmonics of 100 Hz, then undo the resampling
%   to get back to the original pitch.
%   F returns the resampled (flattened) signal, and G is F after
%   filtering.  
%   M is the map from D's timebase to F's.
% 2014-05-01 Dan Ellis dpwe@ee.columbia.edu

doplot = 1;

% (1) Run SAcC pitch tracker with the noisy rats classifier, and
% with the unvoiced state discounted, to get a near-continuous
% pitch track.
[pitch, pvx, times] = sacc_pitchtrack(d, sr);

% Make all frames have a pitch: fill any zeros with .. something
target_pitch = 100.0;
pitch(find(pitch==0)) = target_pitch;

% Build the time mapping
% Map is a two rows; first row indicates a point in original space
% corresponding element in second row is where it ends up
vmap = [times' ; ...
        times(1) + [0, cumsum( pitch(1:end-1)'/target_pitch ...
                               .* diff(times') ) ] ];

% resample waveform to flatten pitches to target_pitch
dm = resample_map(d', sr, vmap);

% Enhance components at the target pitch period
%dmf = enhance_period(dm, round(sr/target_pitch));

%win_t = 15; % median filter time window, in ? 8 ms steps
%win_f = 7;  % local average window in frequency (center pt not used)
%dmf = sgram_enhance(dm, win_t, win_f);

dmf = wiener_icsi(dm, sr);

% Resample back to original domain
du = resample_map(dmf, sr, inv_map(vmap));

% crossfade based on pvx (interpolated up to sr)
do_crossfade = 0;
if do_crossfade
  pvx_dsfact = round(sr*median(diff(times)));
  du = crossfade(du', d, pvx, pvx_dsfact);
end
  
% Return values
y = du;
f = dm;
g = dmf;
m = vmap;

% Plotting carved out to make code more readable

if doplot

  subplot(412)
  sg_fft = 1024;
  sg_olp = sg_fft - sg_fft/8;
  specgram(dm, sg_fft, sr, sg_fft, sg_olp);
  title(['resampled to pitch = ', num2str(target_pitch), ' Hz']);
  cax = [-60 20];
  caxis(cax);
  colormap(1-gray);
  axis([0 length(dm)/sr 0 2000]);

  subplot(413)
  specgram(dmf, sg_fft, sr, sg_fft, sg_olp);
  title('Filtered to enhance 100 Hz');
  caxis(cax);
  colormap(1-gray);
  axis([0 length(dm)/sr 0 2000]);

  subplot(414)
  specgram(du, sg_fft, sr, sg_fft, sg_olp);
  title('Resampled back to original pitch and mixed with original');
  caxis(cax);
  colormap(1-gray);
  axis([0 length(du)/sr 0 2000]);

end

