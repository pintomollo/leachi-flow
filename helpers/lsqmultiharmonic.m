function best = lsqmultiharmonic(x, y, nharm)

  max_iter = 20000;
  damp = 0.2;
  thresh = 0.00000001;
  exp_factor = 3;

  thresh = thresh * damp;
  damp_param = (exp_factor/max_iter);

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
  period_inv = indx(1)/(2*npos);

  goods = (~isnan(x) & ~isnan(y));
  x = x(goods);
  y = y(goods);

  npts = length(x);

  mat = ones(npts, 1+2*nharm);
  freq = 2*pi*x;

  %{
  for i=1:nharm
    curr_freq = freq .* i * period_inv;
    mat(:,2*i) = cos(curr_freq);
    mat(:,2*i+1) = sin(curr_freq);
  end

  params = mat \ y;

  mat = [NaN(npts, 1) mat];
  params = [0; params];
  %}

  for j=1:max_iter
    for i=1:nharm
      curr_freq = freq .* i * period_inv;
      mat(:,2*i) = cos(curr_freq);
      mat(:,2*i+1) = sin(curr_freq);
    end

    params = mat \ y;

    %mat(:,1) = 0;
    counter = zeros(npts, 1);
    for i=1:nharm
      %prev_freq = freq .* i * period_inv;
      %curr_freq = prev_freq;
      %curr_freq = freq .* i * (period_inv + params(1));
      curr_freq = freq .* i * period_inv;
      xi = x * i;

      counter = counter - params(2*i) * xi .* sin(curr_freq) + params(2*i+1) * xi .* cos(curr_freq);
      %mat(:,1) = mat(:,1) - params(2*i+1) * xi .* sin(prev_freq) + params(2*i+2) * xi .* cos(prev_freq);
      %mat(:,2*i+1) = cos(curr_freq);
      %mat(:,2*i+2) = sin(curr_freq);
    end

    params = [counter mat] \ y;

    dw = exp(-j*damp_param)*damp*params(1);
    period_inv = period_inv + dw;

    if (abs(dw) < thresh)
      break;
    end
  end
  %period_inv = period_inv + params(1);

  period = 1/period_inv;
  ampls = sqrt(params(3:2:end).^2 + params(4:2:end).^2);
  phases = atan2(params(4:2:end), params(3:2:end));

  best = [period; ampls; phases];

  return;
end
