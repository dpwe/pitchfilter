function y = enhance_period(x, pd)
% y = enhance_period(x, pd)
%   Filter x to enhance component with period pd.
% 2014-05-14 Dan Ellis dpwe@ee.columbia.edu

% Apply 100 Hz enhancement
%f_ord = round(sr/target_pitch);
f_ord = pd;

if 0
  % Simple comb - feedback at the period, maybe a zero nearby to
  % limit impact
  pole_radius = 0.995;
  zero_radius = 0.1;
  f_a = zeros(1, f_ord+1);
  f_a(1) = 1;
  f_a(f_ord+1) = -(pole_radius^(f_ord));
  f_b(1) = 1;
  f_b(f_ord+1) = -(zero_radius^(f_ord));
  y = filter(f_b, f_a, x);

else
  % Combify an FIR LPF
  fir_ord = 16;
  cutoff = 0.2;
  b = fir1(fir_ord, cutoff);
  % Comb-ize
  f_b = zeros(1, 1 + length(b)*f_ord);
  for i = 1:length(b)
    f_b(1 + (i-1)*f_ord) = b(i);
  end
  % Apply filter
  y = conv(f_b, x);
  % Remove linear phase delay
  y = y(((length(b)-1)/2*f_ord)+[1:length(x)]);
end

% Plot enhancement filter?
%  figure(2)
%  freqz(f_b, f_a);
%  %zplane(f_b, f_a);

