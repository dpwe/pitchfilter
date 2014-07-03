function y = derumble(d, sr, f_hp)
% y = derumble(d, sr, f_hp)
%   Simple HPF removes energy below f_hp (200 Hz by default).
% 2014-07-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;    f_hp = 200;  end

% rumble-removing hp filter is not in icsi's nr, but it's often useful

% Just cut out the low freq for rumble
hp_alpha = 2 * pi * f_hp/sr;
% Zero to remove DC and scale to compensate for pole gain
b_hpf = (1 - hp_alpha/2) * [1 -1];
% Pole to restore low frequencies above f_lp
a_hpf = [1 -(1-hp_alpha)];
% filter
y = filter(b_hpf, a_hpf, d);

%freqz(b_hpf, a_hpf);
  
