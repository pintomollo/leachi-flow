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

    ts = ts/ts(end)

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
  all_data = NaN(2,0);

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

      new_ys = new_ys*sign(acor(I));
    else
      ref = new_ys;
    end

    signals{i,1} = new_ts(goods)+shift;
    signals{i,2} = new_ys;

    plot(signals{i,1}, signals{i,2});

    all_data = [all_data [signals{i,1}; signals{i,2}]];
  end

  all_ts = unique(all_data(1,:));

  %% Need to remove the single data points !!

  [pp1, pthresh] = csaps(all_data(1,:), all_data(2,:));
  %[pp2] = csaps(all_data(1,:), all_data(2,:), pthresh/(dt^3));

  avg_val1 = fnval(pp1, all_ts);
  %avg_val2 = fnval(pp2, all_ts);

  plot(all_ts, avg_val1, 'r')
  %plot(all_ts, avg_val2, 'k')


  keyboard

  return;
end
