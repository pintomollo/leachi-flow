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
  proj_dist = vessel_width / 2;

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

  subplot(2,2,4)
  imagesc(mask);hold on
  plot(branches(:,1), branches(:,2), 'k');

  prev_indx = -1;
  dist = [];

  x1 = branches(1:3:end, 1);
  x2 = branches(2:3:end, 1);
  y1 = branches(1:3:end, 2);
  y2 = branches(2:3:end, 2);

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

  dists(frac < 0 | frac > 1) = Inf;
  inside = (sum(dists < proj_dist, 2)==1);
  mask = reshape(inside, size(X));

  keyboard


  %windows = [32 32; 16 16; 8 8; 8 8];
  %threshs = [Inf; Inf; 10; 5];
  windows = 2.^(max(nextpow2(3*vessel_width)-[0 1 2 2], 2));
  %windows = [64 64; 32 32; 16 16; 16 16];
  threshs = [5; 5; 3; 3];

  data = cell(nframes, 1);
  for nimg=1:nframes-1
    if (prev_indx == nimg)
      img = img_next;
    else
      img = double(load_data(myrecording.channels(1), nimg));
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

    prev_indx = nimg;
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
