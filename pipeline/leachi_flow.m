function [myrecording, opts] =leachi_flow(myrecording, opts)

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

  if (isempty(myrecording.segmentations))
    segmentations = get_struct('segmentation');
  else
    segmentations = myrecording.segmentations;
  end

  b_leachi = get_struct('botrylloides_leachi');

  diff_thresh = 0.75;
  min_length = 1;
  min_branch = 10;
  prop_thresh = 0.3;

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));

  if (isempty(segmentations.detections) || length(segmentations.detections) < nframes)
    detections = get_struct('detection', [1 nframes]);
  else
    detections = segmentations.detections(1:nframes);
  end

  vessel_width = ceil(min(b_leachi.vessel_width.mu) / opts.pixel_size);
  proj_dist = vessel_width * 0.75;

  if (isempty(detections(1).cluster))
    if (opts.verbosity > 1)
      hwait = waitbar(0,'Computing mask...','Name','B. leachi flow');
    end

    orig_img = double(load_data(myrecording.channels(1), 1));

    disk1 = strel('disk', 2*vessel_width);
    disk2 = strel('disk', vessel_width);

    sigma = min(b_leachi.blood_cell.mu(:,1)/(2*opts.pixel_size));

    noise = estimate_noise(orig_img);
    orig_img = gaussian_mex(orig_img, sigma);
    prev_img = padarray(orig_img, [3 3]*vessel_width, NaN);
    %prev_mask = imdilate(imopen(prev_img < noise(1) - diff_thresh*noise(2), disk1), disk2);

    detections(1).noise = noise;

    %figure;
    mask = zeros(img_size+6*vessel_width);
    for nimg=2:nframes
      new_img = double(load_data(myrecording.channels(1), nimg));
      new_img = gaussian_mex(new_img, sigma);

      [img_diff, moire] = immoire(new_img - orig_img, 5, 2*sigma);

      orig_img = new_img;
      img = padarray(orig_img, [3 3]*vessel_width, NaN);
      %curr_mask = imdilate(imopen(img < noise(1) - diff_thresh*noise(2), disk1), disk2);

      img_diff = abs(padarray(img_diff, [3 3]*vessel_width, NaN));
      %img_diff(prev_mask | curr_mask) = false;

      bw = img_diff > diff_thresh * noise(2);
      bw = bwareaopen(bw, ceil(5 / opts.pixel_size).^2);

      %subplot(2,2,1);imagesc(prev_img)
      %subplot(2,2,1);imagesc(prev_mask | curr_mask)
      %subplot(2,2,2);imagesc(img_diff)
      if (any(bw(:)))

        closed = imdilate(bw, disk1);
        open = imerode(closed, disk2);
        mask = mask + open;

        %subplot(2,2,3);imagesc(bw);
        %subplot(2,2,4);imagesc(mask)
        %props = sum(open(:)) * inelems;
        %title(props)
      end

      prev_img = img;
      %prev_mask = curr_mask;

      if (opts.verbosity > 1)
        waitbar(nimg/nframes,hwait);
      end
      %drawnow
    end

    if (opts.verbosity > 1)
      close(hwait);
    end

    keyboard

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

    detections(1).properties = branches;

    %subplot(2,2,4)
    %imagesc(mask);hold on
    %plot(branches(:,1), branches(:,2), 'k');

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

    dists(frac < 0 | frac > 1 | dists > 1) = Inf;
    crosses = any(nodes < 1.5, 2);
    dists(crosses,:) = Inf;
    inside = (sum(isfinite(dists), 2)==1);
    mask = reshape(inside, size(X));
    [junk, indexes] = min(dists(inside,:), [], 2);
    mapping = double(mask);
    mapping(mask) = indexes;

    detections(1).cluster = mapping;
    segmentations.detections = detections;
    myrecording.segmentations = segmentations;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  else
    mapping = detections(1).cluster;
    mask = logical(mapping);

    branches = detections(1).properties;

    %subplot(2,2,4)
    %imagesc(mask);hold on
    %plot(branches(:,1), branches(:,2), 'k');

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
  end

  prev_indx = -1;
  real_mapping = [];

  %windows = [32 32; 16 16; 8 8; 8 8];
  %threshs = [Inf; Inf; 10; 5];
  %windows = 2.^(max(nextpow2(3*vessel_width)-[0 1 2 2], 2));
  windows = vessel_width * [4 3 2 1 1].';
  %windows = [64 64; 32 32; 16 16; 16 16];
  threshs = [5; 5; 3; 3];

  data = cell(nframes-1, 1);

  if (isempty(detections(1).carth) || all(isnan(detections(1).carth(:))))

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
      detections(nimg).carth = speed;

      prev_indx = nimg;
    end

    segmentations.detections = detections;
    myrecording.segmentations = segmentations;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  else
    for nimg=1:nframes-1
      data{nimg} = detections(nimg).carth;
    end
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

  avgs = bsxfun(@times, avgs(:,~empty_branches), sames);
  %{
  figure;
  for i=1:size(avgs,2)
    subplot(1,size(avgs,2), i);
    plot(avgs(:,i));
  end
  %}

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
  %bads = cellfun('isempty', data);
  group_indxs = group_indxs*opts.time_interval;
  [gpos, indxi, indxj] = unique(group_indxs);

  %pos = pos(~bads);
  %pos = pos(:)*opts.time_interval;

  %{
  figure;
  for i=1:size(avgs,2)
    subplot(1,size(avgs,2), i);
    goods = (avgs_indxs == i);

    if (any(goods))
      boxplot(speeds(goods), group_indxs(goods), 'position', gpos);
    end
  end
  %}

  if (opts.verbosity > 1)
    hfig = figure;hold on;
    try
      boxplot(speeds, group_indxs, 'position', gpos);
    catch ME
      keyboard
    end

    navgs = size(avgs, 2);
    colors = redbluemap(navgs);
    for i=1:navgs
      plot(gpos, avgs(:,i), 'Color', colors(i,:), 'LineWidth', 2);
    end
  end

  goods = (~isnan(group_indxs) & ~isnan(speeds));
  prev_params = -Inf;
  for i=1:20
    [period, ampls, phases] = lsqmultiharmonic(group_indxs(goods), speeds(goods));
    nharm = length(ampls);

    sign_val = zeros(size(gpos));
    for j=1:nharm
      sign_val = sign_val + ampls(j)*cos((j*gpos/period)*2*pi + phases(j));
    end
    thresh = max(ampls)/2;

    if (opts.verbosity > 1)
      plot(gpos, sign_val, 'k');
    end

    goods = (speeds < sign_val(indxj) + thresh & speeds > sign_val(indxj) - thresh);

    dx = sum([period ampls(1) phases(1)] - prev_params);
    prev_params = [period ampls(1) phases(1)];

    if (dx < 1e-6)
      break;
    end
  end

  if (isempty(myrecording.trackings))
    trackings = get_struct('tracking');
  else
    trackings = myrecording.trackings;
  end

  if (isempty(trackings.detections))
    res = get_struct('detection', 1);
  else
    res = trackings.detections(1);
  end

  res.carth = [speeds group_indxs];
  res.cluster = goods;
  res.properties = [ampls(:).'; period*ones(1,length(ampls)); phases(:).'];

  if (opts.verbosity > 1)
    plot(gpos, sign_val, 'k', 'LineWidth', 2);
    tmp = 2*sum(ampls(:));
    ylim([-tmp tmp]);
    title(num2str(res.properties));
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_f.png'], '-dpng');
    delete(hfig);
  end

  trackings.detections = res;

  myrecording.segmentations = segmentations;
  myrecording.trackings = trackings;

  return;
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
