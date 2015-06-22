function lsqmultiharmonic(x, y, nharm)

  x = x(:);
  y = y(:);

  pos = unique(x);
  dt = min(diff(pos));

  pos = pos / dt;
  x = x / dt;

  uniform_pos = [min(pos):max(pos)];
  npos = length(uniform_pos);

  mean_val = mymean(y(:), 1, x(:));

  goods = (~isnan(pos) & ~isnan(mean_val));

  mean_val = interp1(pos(goods), mean_val(goods), uniform_pos);

  nfft = 2^nextpow2(npos);

  ff_val = fft(mean_val, nfft)/npos;
  ff_val = 2*abs(ff_val(1:nfft/2+1));

  [junk, indx] = max(ff_val);
  period = 2*npos/indx(1);

  goods = (~isnan(x) & ~isnan(y));
  x = x(goods);
  y = y(goods);

  npts = length(x);

  mat = ones(npts, 1+2*nharm);
  freq = 2*pi*x;

  for i=1:nharm
    curr_freq = freq .* i/period;
    mat(:,2*i) = cos(curr_freq);
    mat(:,2*i+1) = sin(curr_freq);
  end

  params = mat \ y;

  keyboard

  return;
end
