function flows = combine_flows(flows)

  if (ischar(flows))
    names = dir(flows);
    nflows = length(names);

    flows = cell(nflows, 2);
    for i=1:nflows
      data = load(names(i).name);
      flow = data.myrecording.trackings.detections.properties;
      time = [1:size_data(data.myrecording.channels(1))] * data.opts.time_interval;

      flows{i,1} = flow;
      flows{i,2} = time;
    end
  end

  nflows = size(flows, 1);
  figure;hold on;

  signals = cell(nflows, 2);
  %all_t = NaN(1,0);
  dt = Inf;
  tmin = Inf;
  tmax = -Inf;

  for i=1:nflows
    ts = flows{i,2};
    ys = fnval(flows{i,1}, ts);

    %ts = ts/ts(end)

    dt = min(dt, median(diff(ts)));
    tmin = min(tmin, min(ts));
    tmax = max(tmax, max(ts));

    plot(ts, ys);

    %all_t = [all_t ts];
    signals{i,1} = ts;
    signals{i,2} = ys;
  end

  %all_t = unique(all_t);
  %dt = median(diff(all_t));

  new_ts = [tmin:dt:tmax];

  params = realign_curves(new_ts, signals, []);
  %{
  all_data = NaN(1,0);

  params = zeros(nflows, 2);

  figure;hold on;
  for i=1:nflows
    ts = signals{i,1};
    ys = signals{i,2};

    goods = (new_ts >= min(ts) & new_ts <= max(ts));
    new_ys = interp1(ts, ys, new_ts(goods));
    shift = 0;

    if (i > 1)
      [acor,lag] = xcorr(ref,new_ys);

      [~,I] = max(abs(acor));
      lagDiff = lag(I);
      shift = lagDiff*dt;

      mysign = sign(acor(I));
    else
      ref = new_ys;
      mysign = 1;
    end

%    signals{i,1} = new_ts(goods)+shift;
%    signals{i,2} = new_ys;

%    plot(signals{i,1}, signals{i,2});
    plot(new_ts(goods)+shift, mysign*new_ys);
    params(i,:) = [shift mysign];

    all_data = [all_data new_ts(goods)+shift];
  end
  %}

  new_ts = [params(end,1):dt:params(end,2)];
  all_data = NaN(nflows, length(new_ts));

  figure;hold on;
  for i=1:nflows
    ts = signals{i,1};
    f = flows{i,1};

    s = params(i,1);
    goods = (new_ts-s >= min(ts) & new_ts-s <= max(ts));
    ys = params(i,2)*fnval(f, new_ts(goods) - s);

    all_data(i,goods) = ys;

    plot(new_ts, all_data(i,:))
  end

  %% Need to remove the single data points !!

  %[pp1, pthresh] = csaps(all_data(1,:), all_data(2,:));
  %[pp2] = csaps(all_data(1,:), all_data(2,:), pthresh/(dt^3));

  %avg_val1 = fnval(pp1, all_ts);
  %avg_val2 = fnval(pp2, all_ts);

  %plot(all_ts, avg_val1, 'r')
  %plot(all_ts, avg_val2, 'k')
  [m,s,n] = mymean(all_data);
  goods = (n>2);

  all_ts = repmat(new_ts(goods), 1, nflows);
  all_ys = all_data(:, goods).';

  [pp] = csaps(all_ts(:), all_ys(:), 1/max(all_ts));
  new_ts = new_ts(goods);
  ref = fnval(pp, new_ts);

  params2 = realign_curves(new_ts, signals, ref);

  new_ts = [params2(end,1):dt:params2(end,2)];
  all_data = NaN(nflows, length(new_ts));

  figure;hold on;
  for i=1:nflows
    ts = signals{i,1};
    f = flows{i,1};

    s = params2(i,1);
    goods = (new_ts-s >= min(ts) & new_ts-s <= max(ts));
    ys = params2(i,2)*fnval(f, new_ts(goods) - s);

    all_data(i,goods) = ys;

    plot(new_ts, all_data(i,:))
  end

  [m,s,n] = mymean(all_data);
  goods = (n>3);

  all_ts = repmat(new_ts(goods), 1, nflows);
  all_ys = all_data(:, goods).';

  [pp] = csaps(all_ts(:), all_ys(:), 1/max(all_ts));
  new_ts = new_ts(goods);
  ref2 = fnval(pp, new_ts);

  cross = find(ref2(2:end) > 0 & ref2(1:end-1) <= 0);

  if (length(cross) > 1)

    overlap = min(200, min(cross(1)-1, size(all_data,2) - cross(2)));

    before = all_data(:,1:cross(1)-1+overlap);
    data = all_data(:,cross(1)-overlap:cross(2)-1+overlap);
    after = all_data(:,cross(2)-overlap:end);

    npts = size(data,2);
    all_data = NaN(3*nflows, npts);
    all_data(1:nflows, end-size(before,2)+1:end) = before;
    all_data([1:nflows]+nflows, :) = data;
    all_data([1:nflows]+2*nflows, 1:min(size(after,2),end)) = after(:,1:min(end,npts));

    orig = new_ts(cross(1));
    new_ts = new_ts(cross(1)-overlap:cross(2)-1+overlap);
  elseif (length(cross) == 1)
    orig = new_ts(cross(1));
  else
    orig = 0;
  end


  [m,s,n] = mymean(all_data);
  goods = (n>3);

  keyboard

  all_ts = repmat(new_ts(goods), 3*nflows, 1).';
  all_ys = all_data(:, goods).';

  [pp] = csaps(all_ts(:), all_ys(:), 1/max(all_ts(:)));
  ref3 = fnval(pp, new_ts(goods));
  stds = sqrt(mymean(bsxfun(@minus, all_ys.', ref3).^2, 1));

  figure;plot(repmat(new_ts - orig, 3*nflows, 1).', all_data.', 'b');
  hold on;
  plot(new_ts(goods) - orig, ref3, 'r')
  plot(new_ts(goods) - orig, ref3+stds, 'k')
  plot(new_ts(goods) - orig, ref3-stds, 'k')

  return;
end

function [params] = realign_curves(new_ts, signals, ref)

  nflows = size(signals, 1);
  dt = median(diff(new_ts));

  all_data = NaN(1,0);
  params = zeros(nflows+1, 2);

  figure;hold on;
  for i=1:nflows
    ts = signals{i,1};
    ys = signals{i,2};

    goods = (new_ts >= min(ts) & new_ts <= max(ts));
    new_ys = interp1(ts, ys, new_ts(goods));
    shift = 0;

    if (isempty(ref))
      ref = new_ys;
      mysign = 1;
    else
      [acor,lag] = xcorr(ref,new_ys);

      [~,I] = max(abs(acor));
      lagDiff = lag(I);
      shift = lagDiff*dt;

      mysign = sign(acor(I));
    end

%    signals{i,1} = new_ts(goods)+shift;
%    signals{i,2} = new_ys;

%    plot(signals{i,1}, signals{i,2});
    plot(new_ts(goods)+shift, mysign*new_ys);
    params(i,:) = [shift mysign];

    all_data = [all_data new_ts(goods)+shift];
  end

  params(end,:) = [min(all_data) max(all_data)];

  return;
end
