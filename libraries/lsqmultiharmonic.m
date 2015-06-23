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

  %%%% interpolated FFT from ACD test toolbox (Tamas Virostek)
  Fc=fft(mean_val);
  F=abs(Fc);
  [Mfft,w]=max(F(2:round(npos/2)));

  if w>1
      %calculating the 2 points, between them the estimated frequency is
      if F(w-1)>F(w+1) w=w-1; end
  end

  n=2*pi/npos;
  U=real(Fc(w+1));    V=imag(Fc(w+1));
  U1=real(Fc(w+2));  V1=imag(Fc(w+2));
  Kopt=(sin(n*w)*(V1-V)+cos(n*w)*(U1-U))/(U1-U);
  Z1=V*(Kopt-cos(n*w))/sin(n*w)+U;
  Z2=V1*(Kopt-cos(n*(w+1)))/sin(n*(w+1))+U1;

  lambda=acos((Z2*cos(n*(w+1))-Z1*cos(n*w))/(Z2-Z1))/n;
  period_inv=lambda/npos;
  %%%%

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
  phases = atan2(-params(4:2:end), params(3:2:end));

  best = [period; ampls; phases];

  return;
end
