function leachi_flow(myrecording, opts)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

  if (~isstruct(myrecording))
    [myrecording, opts] = inspect_recording();

    if (~isempty(myrecording))
      [myrecording, opts] = preprocess_movie(myrecording, opts);
      save(myrecording.experiment, 'myrecording', 'opts');
    end
  end

  diff_thresh = 5;
  min_length = 1;
  min_branch = 10;
  prop_thresh = 0.5;

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));

  orig_img = double(load_data(myrecording.channels(1), 1));
  %prev_noise = estimate_noise(prev_img);

  vessel_width = ceil(15 / opts.pixel_size);
  disk1 = strel('disk', 2*vessel_width);
  disk2 = strel('disk', vessel_width);

  noise = estimate_noise(orig_img);
  prev_img = padarray(orig_img, [3 3]*vessel_width);
  prev_mask = imdilate(imopen(prev_img < noise(1) + diff_thresh*noise(2), disk1), disk2);

  indexes = NaN(nframes, 2);
  masks = cell(nframes, 1);
  indx = 1;

  %figure;
  mask = zeros(img_size+6*vessel_width);
  for nimg=2:nframes
    orig_img = double(load_data(myrecording.channels(1), nimg));
    img = padarray(orig_img, [3 3]*vessel_width);
    curr_mask = imdilate(imopen(img < noise(1) + diff_thresh*noise(2), disk1), disk2);

    %bkg_diff = abs(prev_noise(1) - noise(1));
    %curr_noise = max(prev_noise(2), noise(2));

    img_diff = abs(prev_img - img);

    img_diff(prev_mask | curr_mask) = false;

    %bw = img_diff > bkg_diff + diff_thresh * curr_noise;
    bw = img_diff > diff_thresh * noise(2);

    %if (nimg==5)
    %figure;imagesc(bw)
    %end

    bw = bwareaopen(bw, ceil(5 / opts.pixel_size).^2);

    if (any(bw(:)))

      closed = imdilate(bw, disk1);
      open = imerode(closed, disk2);

      props = sum(open(:)) * inelems;

      if (props > 0.5)
        masks{indx} = mask;
        mask = zeros(img_size+6*vessel_width);
        indx = nimg - 1;
      else
        mask = mask + open;
      end
      indexes(nimg - 1, 1) = indx;
      %closed = imclose(bw, disk);
      %open = imopen(closed, disk);

      %subplot(2,2,1);imagesc(prev_img)
      %subplot(2,2,2);imagesc(closed)
      %subplot(2,2,3);imagesc(open);
      %subplot(2,2,4);imagesc(mask)
      %title(props)

    else
      masks{indx} = mask;
      mask = zeros(img_size+6*vessel_width);

      noise = estimate_noise(orig_img);
      indx = nimg - 1;
    end

    prev_img = img;
    %prev_noise = noise;
    prev_mask = curr_mask;

    
    %drawnow
  end
  masks{indx} = mask;

  %keyboard

  [sorted, coords, inverse] = unique(indexes(:,1));
  [tmp_coords, tmp_indx] = sort(coords);
  lens = diff([tmp_coords; nframes]);
  [junk,rev_indx] = sort(tmp_indx);
  goods = (lens(rev_indx)>min_length) & isfinite(sorted);
  groups = sorted(goods);

  ngroups = length(groups);

  indexes(~ismember(indexes, groups), 1) = NaN;
  indexes(groups, 2) = [1:ngroups];
  masks = masks(groups);
  centers = cell(size(masks));
  %proj_dist = 8;
  proj_dist = vessel_width / 2;

  %keyboard

  for i=1:ngroups
    mask = masks{i};
    %mask = mask(3*vessel_width+[1:img_size(1)], 3*vessel_width+[1:img_size(2)]);
    mask = imnorm(mask);

    figure;subplot(2,2,1)
    imagesc(mask)

    mask = (mask > prop_thresh);

    subplot(2,2,2)
    imagesc(mask)

    mask = bwareaopen(mask, vessel_width^2);

    subplot(2,2,3)
    imagesc(mask)

    mask = bwmorph(mask, 'thin', Inf);
    [icoord, jcoord] = find(mask);
    centers{i} = [jcoord, icoord];
    tmp_mask = imdilate(mask, disk2);
    masks{i} = tmp_mask(3*vessel_width+[1:img_size(1)], 3*vessel_width+[1:img_size(2)]);

    subplot(2,2,4)
    imagesc(tmp_mask+mask)
  end
  prev_indx = -1;

  keyboard

  %windows = [32 32; 16 16; 8 8; 8 8];
  %threshs = [Inf; Inf; 10; 5];
  windows = 2.^(max(nextpow2(3*vessel_width)-[0 1 2 2], 2));
  %windows = [64 64; 32 32; 16 16; 16 16];
  threshs = [5; 5; 3; 3];

  data = cell(nframes, 1);

  for nimg=1:nframes-1
    if (~isnan(indexes(nimg, 1)))
      if (prev_indx == nimg)
        img = img_next;
      else
        img = double(load_data(myrecording.channels(1), nimg));
        mask = masks{indexes(indexes(nimg, 1), 2)};
        dist = [];

        branches = sort_shape(centers{indexes(indexes(nimg, 1), 2)});
        branches_size = cellfun(@(x)(size(x,1)), branches);
        goods = (branches_size > min_branch);

        if (~any(goods))
          continue;
        end

        branches = branches(goods);
        branches_size = branches_size(goods);
        center = cat(1, branches{:});

        vectors = cellfun(@(x)(differentiator(carth2linear(x), x, 'replicate', 1, 11)), branches, 'UniformOutput', false);
        parallels = cat(1,vectors{:});
        parallels = parallels(:,1:2);
        parallels = bsxfun(@rdivide, parallels, sqrt(sum(parallels.^2, 2)));
      end
      img_next = double(load_data(myrecording.channels(1), nimg+1));

      [x,y,u,v] = matpiv_nfft(img, img_next, windows, 1/32, threshs, mask);
      empties = (u == 0 & v == 0);
      u(empties) = NaN;
      v(empties) = NaN;

      if (isempty(dist) && ~isempty(center))
        dist = sqrt(bsxfun(@minus, center(:,1), x(:).').^2 + ...
                    bsxfun(@minus, center(:,2), y(:).').^2);
        dist(dist>3*proj_dist) = Inf;

        weights = exp(-(dist.^2)/(2*((proj_dist/2)^2)));

        weights = bsxfun(@rdivide, weights, sum(weights, 2));
      end

      %subplot(2,2,1);imagesc(img)
      %subplot(2,2,2);imagesc(img_next)
      %subplot(2,2,3);imagesc(mask);
      %subplot(2,2,4);quiver(x,-y,u,-v)
      %axis([1 img_size(2) -img_size(1) -1])

      if (~isempty(weights))
        speed_x = nansum(bsxfun(@times, weights, u(:).'), 2);
        speed_y = nansum(bsxfun(@times, weights, v(:).'), 2);

        speed = dot(parallels, [speed_x speed_y], 2);

        results = mat2cell(speed, branches_size, 1);
        data{nimg} = results;
      end
      %keyboard

      prev_indx = nimg;
    end
  end

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

  speeds = [];
  group_indxs = [];
  pos = [1:nframes];
  for i = pos
    tmp = data{i};
    speeds = [speeds; tmp];
    group_indxs = [group_indxs; ones(size(tmp))*i];
  end
  bads = cellfun('isempty', data);
  pos = pos(~bads);

  figure;boxplot(speeds, group_indxs, 'position', pos);

  return;
end
