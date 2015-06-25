function [period, ampls, phases] = lsqmultiharmonic(x, y, nharm)

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

  if (nargin < 3)
    nharm = estimate_harmonics_number(mean_val);
    nharm = max(nharm, 1);
  end

  %%%% interpolated FFT from ACD test toolbox (Tamas Virostek)
  Fc=fft(mean_val);
  F=abs(Fc(2:round(npos/2)));

  %figure;plot(F);hold on;
  if (nharm < 2)
    [Mfft,w]=max(F);
  else

    [maxs, indxs] = local_extrema(F);
    [maxs, s] = sort(maxs, 'descend');
    indxs = indxs(s);

    %scatter(indxs, maxs, 'g')

    if (nharm == 2)
      w = min(indxs(1:2));
    else
      dist = abs(indxs(2:end) - indxs(1));
      ratio = dist(1) ./ dist;
      new_val = round(ratio)*dist(1);

      ratio_inv = dist./dist(1);
      new_val2 = round(ratio).*dist;

      new_val(ratio < 1) = new_val2(ratio < 1);

      goods = [true (abs(dist - new_val) < 2)];

      w = min(indxs(goods));

      %scatter(indxs(goods), F(indxs(goods)), 'r')
    end
  end
  %scatter(w, F(w), 'k')

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

  for j=1:max_iter
    for i=1:nharm
      curr_freq = freq .* i * period_inv;
      mat(:,2*i) = cos(curr_freq);
      mat(:,2*i+1) = sin(curr_freq);
    end

    params = mat \ y;

    counter = zeros(npts, 1);
    for i=1:nharm
      curr_freq = freq .* i * period_inv;
      xi = x * i;

      counter = counter - params(2*i) * xi .* sin(curr_freq) + params(2*i+1) * xi .* cos(curr_freq);
    end

    params = [counter mat] \ y;

    dw = exp(-j*damp_param)*damp*params(1);
    period_inv = period_inv + dw;

    if (abs(dw) < thresh)
      break;
    end
  end

  period = 1/period_inv;
  ampls = sqrt(params(3:2:end).^2 + params(4:2:end).^2);
  phases = atan2(-params(4:2:end), params(3:2:end));

  return;
end
