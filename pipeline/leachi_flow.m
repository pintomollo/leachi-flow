function leachi_flow(myrecording, opts)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

  if (~isstruct(myrecording))
    if (ischar(myrecording))
      [myrecording, opts] = inspect_recording(myrecording);
    else
      [myrecording, opts] = inspect_recording();
    end

    if (~isempty(myrecording))
      [myrecording, opts] = preprocess_movie(myrecording, opts);
      save([myrecording.experiment '.mat'], 'myrecording', 'opts');
    end
  end

  b_leachi = get_struct('botrylloides_leachi');

  diff_thresh = 2;
  min_length = 1;
  min_branch = 10;
  prop_thresh = 0.3;

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));

  orig_img = double(load_data(myrecording.channels(1), 1));

  vessel_width = ceil(min(b_leachi.vessel_width.mu) / opts.pixel_size);
  proj_dist = vessel_width * 0.75;

  disk1 = strel('disk', 2*vessel_width);
  disk2 = strel('disk', vessel_width);

  sigma = min(b_leachi.blood_cell.mu(:,1)/(2*opts.pixel_size));

  noise = estimate_noise(orig_img);
  prev_img = padarray(gaussian_mex(orig_img, sigma), [3 3]*vessel_width, NaN);
  prev_mask = imdilate(imopen(prev_img < noise(1) - diff_thresh*noise(2), disk1), disk2);

  %figure;
  mask = zeros(img_size+6*vessel_width);
  for nimg=2:nframes
    orig_img = double(load_data(myrecording.channels(1), nimg));
    img = padarray(gaussian_mex(orig_img, sigma), [3 3]*vessel_width, NaN);
    curr_mask = imdilate(imopen(img < noise(1) - diff_thresh*noise(2), disk1), disk2);

    img_diff = abs(prev_img - img);
    img_diff(prev_mask | curr_mask) = false;

    bw = img_diff > diff_thresh * noise(2);
    bw = bwareaopen(bw, ceil(5 / opts.pixel_size).^2);

    if (any(bw(:)))

      closed = imdilate(bw, disk1);
      open = imerode(closed, disk2);
      mask = mask + open;

      %subplot(2,2,1);imagesc(prev_img)
      %subplot(2,2,2);imagesc(closed)
      %subplot(2,2,3);imagesc(open);
      %subplot(2,2,4);imagesc(mask)
      %props = sum(open(:)) * inelems;
      %title(props)
    end

    prev_img = img;
    prev_mask = curr_mask;
    %drawnow
  end

  mask = imnorm(mask);

  %figure;subplot(2,2,1)
  %imagesc(mask)

  mask = (mask > prop_thresh);

  %subplot(2,2,2)
  %imagesc(mask)

  mask = bwareaopen(mask, vessel_width^2);

  %subplot(2,2,3)
  %imagesc(mask)

  mask = bwmorph(mask, 'thin', Inf);
  mask = mask(3*vessel_width+[1:img_size(1)], 3*vessel_width+[1:img_size(2)]);
  [icoord, jcoord] = find(mask);
  centers = [jcoord, icoord];

  branches = sort_shape(centers, min_branch);
  if (isempty(branches))
    error('nothing');
  end

  %subplot(2,2,4)
  %imagesc(mask);hold on
  %plot(branches(:,1), branches(:,2), 'k');

  prev_indx = -1;
  dist = [];

  x1 = branches(1:3:end, 1);
  x2 = branches(2:3:end, 1);
  y1 = branches(1:3:end, 2);
  y2 = branches(2:3:end, 2);

  branches_size = length(x1);

  vects = [x2-x1 y2-y1];
  lens = 1 ./ sum(vects.^2, 2);
  cross = (x2.*y1 - y2.*x1);
  widths = 1/(proj_dist)^2;

  origin = [x1 y1].';
  params = [vects cross lens sqrt(lens .* widths)].';

  [X,Y] = meshgrid([1:img_size(2)], [1:img_size(1)]);

  dists = bsxfun(@times, (bsxfun(@plus, ...
                            bsxfun(@times, X(:), params(2, :)) - ...
                              bsxfun(@times, Y(:), params(1,:)), ...
                            params(3,:))).^2, ...
                         params(4,:) .* widths);

  nodes = widths*(bsxfun(@minus, X(:), [x1;x2].').^2 + ...
                  bsxfun(@minus, Y(:), [y1;y2].').^2);

  frac = bsxfun(@times, bsxfun(@minus, X(:), ...
                                       origin(1,:)), ...
                        params(1,:) .* params(4,:)) + ...
         bsxfun(@times, bsxfun(@minus, Y(:), ...
                                       origin(2,:)), ...
                        params(2,:) .* params(4,:));

  %{
  perps = bsxfun(@times, bsxfun(@minus, X(:), ...
                                       origin(1,:)), ...
                        -params(2,:) .* params(5,:)) + ...
         bsxfun(@times, bsxfun(@minus, Y(:), ...
                                       origin(2,:)), ...
                        params(1,:) .* params(5,:));
  %}

  dists(frac < 0 | frac > 1 | dists > 1) = Inf;
  crosses = any(nodes < 1.5, 2);
  dists(crosses,:) = Inf;
  inside = (sum(isfinite(dists), 2)==1);
  mask = reshape(inside, size(X));
  [junk, indexes] = min(dists(inside,:), [], 2);
  mapping = double(mask);
  mapping(mask) = indexes;

  real_mapping = [];

  %windows = [32 32; 16 16; 8 8; 8 8];
  %threshs = [Inf; Inf; 10; 5];
  %windows = 2.^(max(nextpow2(3*vessel_width)-[0 1 2 2], 2));
  windows = vessel_width * [4 3 2 1 1].';
  %windows = [64 64; 32 32; 16 16; 16 16];
  threshs = [5; 5; 3; 3];

  data = cell(nframes-1, 1);
  for nimg=1:nframes-1
    if (prev_indx == nimg)
      img = img_next;
    else
      img = double(load_data(myrecording.channels(1), nimg));
    end
    img_next = double(load_data(myrecording.channels(1), nimg+1));

    %%%%%%% COULD FILTER OUT VECTORS THAT ARE NOT // WITH THE CENTERS. EITHER DURING OR AFTER THE PIV

    [x,y,u,v] = matpiv_nfft(img, img_next, windows, 1/32, threshs, mask);
    empties = (u == 0 & v == 0);
    u(empties) = NaN;
    v(empties) = NaN;

    if (isempty(real_mapping))
      tmp_vals = bilinear_mex(mapping, x, y);
      real_mapping = blockproc(tmp_vals, [1 1], @local_mapping, 'BorderSize', [1 1], 'TrimBorder', false);

      inside = (real_mapping > 0);

      ngoods = sum(inside(:));
      subs = sub2ind([ngoods branches_size], [1:ngoods].', real_mapping(inside));

      others = ~ismember([1:ngoods*branches_size], subs);
      one = ones(1, branches_size);
    end

    %subplot(2,2,1);imagesc(img)
    %subplot(2,2,2);imagesc(img_next)
    %subplot(2,2,3);imagesc(mask);
    %subplot(2,2,4);quiver(x,-y,u,-v)
    %axis([1 img_size(2) -img_size(1) -1])
    %drawnow

    speed = bsxfun(@times, u(inside), params(1,:) .* sqrt(params(4,:))) + ...
             bsxfun(@times, v(inside), params(2,:) .* sqrt(params(4,:)));

    speed(others) = NaN;

    data{nimg} = speed;

    prev_indx = nimg;
  end

  ndata = length(data);
  avgs = cellfun(@nanmedian, data, 'UniformOutput', false);
  avgs = cat(1, avgs{:});
  corrs = corr(avgs);

  C = corrs - eye(branches_size);
  sC = (C >= 0);
  aC = abs(C);
  groups = NaN(branches_size);
  groups(:,1) = [1:branches_size];

  empty_branches = all(isnan(C), 1);

  for i=1:branches_size^2
    [val, indxi] = max(aC,[],1);
    [val, indxj] = max(val,[],2);

    if (val == 0)
      break;
    end

    indxi = indxi(indxj);

    aC(indxi, indxj) = 0;
    aC(indxj, indxi) = 0;

    rowi = any(abs(groups)==indxi,2);
    rowj = any(abs(groups)==indxj,2);

    if (~any(rowj & rowi))
      posi = any(groups(rowi,:)==indxi);
      posj = any(groups(rowj,:)==indxj);
      posc = sC(indxi, indxj);
      fact = (-1)^(posi+posj+posc+1);

      first = find(isnan(groups(rowi,:)), 1, 'first');
      last = find(~isnan(groups(rowj,:)), 1, 'last');

      groups(rowi, [first:first+last-1]) = fact*groups(rowj, [1:last]);
      groups(rowj,:) = NaN;
    end

    if (sum(~isnan(groups(:,1)))==1)
      break;
    end
  end

  values = groups(~isnan(groups));
  [junk, indx] = sort(abs(values));
  sames = sign(values(indx));
  sames = sames(~empty_branches);
  sames = sames(:).';

  %{
  keyboard

  for i=1:ngroups
    results = data(indexes(1:nframes,1)==groups(i));
    ndata = length(results);
    nsubs = length(results{1});

    %% Might want to use some statistics to determine relevant averages
    %signs = NaN(ndata, nsubs, 4);
    signs = NaN(ndata, nsubs);

    for j=1:ndata
      %sm1 = cellfun(@mean, results{j});
      %sm2 = cellfun(@median, results{j});
      %sm3 = cellfun(@std, results{j});
      %sm4 = cellfun(@(x)(1.4826*mad(x,1)), results{j});

      %signs(j, :, 1) = sm1;
      %signs(j, :, 2) = sm2;
      %signs(j, :, 3) = sm3;
      %signs(j, :, 4) = sm4;

      signs(j,:) = cellfun(@mean, results{j});
    end

    temp_var = sign(mean(signs, 2));
    invert = sign(mean(bsxfun(@times, signs, temp_var)));

    for j=1:ndata
      for k=1:nsubs
        results{j}{k} = results{j}{k} * invert(k);
      end
      results{j} = cat(1, results{j}{:});
    end

    data(indexes(1:nframes,1)==groups(i)) = results;
  end
  %}

  avgs = bsxfun(@times, avgs(:,~empty_branches), sames);
  figure;
  for i=1:size(avgs,2)
    subplot(1,size(avgs,2), i);
    plot(avgs(:,i));
  end

  speeds = [];
  group_indxs = [];
  avgs_indxs = [];
  pos = [1:ndata];
  for i = pos
    tmp_all = bsxfun(@times, data{i}(:,~empty_branches), sames);
    for j=1:size(avgs,2)
      tmp = tmp_all(:,j);
      tmp = tmp(isfinite(tmp));
      speeds = [speeds; tmp];
      group_indxs = [group_indxs; ones(size(tmp))*i];
      avgs_indxs = [avgs_indxs; ones(size(tmp))*j];
    end
  end
  bads = cellfun('isempty', data);
  pos = pos(~bads);
  pos = pos(:)*opts.time_interval;

  figure;
  for i=1:size(avgs,2)
    subplot(1,size(avgs,2), i);
    goods = (avgs_indxs == i);

    if (any(goods))
      boxplot(speeds(goods), group_indxs(goods), 'position', pos);
    end
  end

  figure;boxplot(speeds, group_indxs, 'position', pos);

  [gpos, indxi, indxj] = unique(group_indxs);

  goods = (~isnan(group_indxs) & ~isnan(speeds));
  prev_params = -Inf;
  for i=1:10
    vals = lsqmultiharmonic(group_indxs(goods), speeds(goods), 1);
    bparams = vals([2 1 (end-1)/2+2]);
    hold on;plot(pos, bparams(1)*cos(((pos/bparams(2))*2*pi + bparams(3))), 'k');

    sign_val = bparams(1)*cos(((gpos/bparams(2))*2*pi + bparams(3)));
    thresh = bparams(1)/2;

    goods = (speeds < sign_val(indxj) + thresh & speeds > sign_val(indxj) - thresh);

    dx = sum(bparams - prev_params);
    prev_params = bparams;

    if (dx < 1e-6)
      break;
    end
  end

  keyboard

  avg_avg = mymean(speeds, 1, group_indxs);
  ampl = sqrt(2)*std(speeds);
  freq = 2*pi / (sqrt(2)*std(differentiator(pos, avg_avg/ampl)));
  smooth = differentiator(pos,cumsum(avg_avg));
  cross = find(smooth(1:end-1)>0 & smooth(2:end)<=0);
  if (isempty(cross))
    step = 0;
  else
    step = mean(mod(cross/freq, 1)) - 0.5;

    if (step < 0)
      step = 1 + step;
    end
  end

  p0 = [ampl, freq, step];
  x = group_indxs*opts.time_interval;
  y = speeds;
  [ym, ys] = mymean(y);
  w = exp(-(y - (ym + ys)).^2 / (ys^2)) + exp(-(y - (ym - ys)).^2 / (ys^2));

  [b,f] = myfit(@err, p0, [0 Inf; 1 Inf; 0 1]);

%%  fit_opt = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
%%  best = lsqcurvefit(@sinusoidal_fit, init, group_indxs*opts.time_interval, speeds, [0 opts.time_interval 0], [Inf Inf 1], fit_opt);

  keyboard

  return;

  function val = err(p, junk)

    val = sum((y - p(1)*sin(((x/p(2)) - p(3))*2*pi) .* w).^2);
    hold off;
    scatter(x,y,'r');
    hold on;
    scatter(x,p(1)*sin(((x/p(2)) - p(3))*2*pi),'b');

    drawnow

    return;
  end
end

function index = local_mapping(block)

  index = NaN;

  data = block.data(:);
  center = data(5);
  data = data(data~=0 & mod(data,1)==0);

  if (~isempty(data) && center~=0)
    vals = unique(data);
    if (numel(vals)==1)
      index = vals;
    else
      counts = sum(bsxfun(@eq, vals(:), data(:).'), 2);
      [indxs] = find(counts==max(counts));

      if (length(indxs)==1)
        index = vals(indxs);
      end
    end
  end

  return;
end

function y = sinusoidal_fit(params, x)

  y = params(1)*sin(((x/params(2)) - params(3))*2*pi);

  return;
end
