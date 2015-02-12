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
  prop_thresh = 0.75;

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));

  prev_img = double(load_data(myrecording.channels(1), 1));
  prev_noise = estimate_noise(prev_img);

  disk1 = strel('disk', 40);
  disk2 = strel('disk', 20);

  indexes = NaN(nframes, 2);
  masks = cell(nframes, 1);
  indx = 1;

  %figure;
  mask = zeros(img_size);
  noise = estimate_noise(prev_img);
  for nimg=2:nframes
    img = double(load_data(myrecording.channels(1), nimg));

    %bkg_diff = abs(prev_noise(1) - noise(1));
    %curr_noise = max(prev_noise(2), noise(2));

    img_diff = abs(prev_img - img);
    %bw = img_diff > bkg_diff + diff_thresh * curr_noise;
    bw = img_diff > diff_thresh * noise(2);
    bw = bwareaopen(bw, 4);

    if (any(bw(:)))

      closed = imdilate(bw, disk1);
      open = imerode(closed, disk2);

      props = sum(open(:)) * inelems;

      if (props > 0.25)
        masks{indx} = mask;
        mask = zeros(img_size);
        indx = nimg - 1;
      else
        mask = mask + open;
      end
      indexes(nimg - 1, 1) = indx;
      %closed = imclose(bw, disk);
      %open = imopen(closed, disk);

    %  subplot(2,2,1);imagesc(prev_img)
    %  subplot(2,2,2);imagesc(closed)
    %  subplot(2,2,3);imagesc(open);
    %  subplot(2,2,4);imagesc(mask)
    %  title(props)

    else
      noise = estimate_noise(prev_img);
      indx = nimg - 1;
    end

    prev_img = img;
    prev_noise = noise;

    %drawnow
    %keyboard
  end

  [sorted, coords, inverse] = unique(indexes(:,1));
  goods = [diff(coords)>min_length; false] & isfinite(sorted);
  groups = sorted(goods);

  ngroups = length(groups);

  indexes(~ismember(indexes, groups), 1) = NaN;
  indexes(groups, 2) = [1:ngroups];
  masks = masks(groups);

  for i=1:ngroups
    mask = imnorm(masks{i});
    mask = (mask > prop_thresh);
    mask = bwareaopen(mask, 300);
    mask = bwmorph(mask, 'thin', Inf);
    masks{i} = imdilate(mask, disk2);
  end

  figure;
  prev_indx = -1;

  windows = [32 32; 16 16; 8 8; 8 8];
  threshs = [Inf; Inf; 10; 5];

  for nimg=1:nframes-1
    if (~isnan(indexes(nimg, 1)))
      if (prev_indx == nimg)
        img = img_next;
      else
        img = double(load_data(myrecording.channels(1), nimg));
        mask = masks{indexes(indexes(nimg, 1), 2)};
      end
      img_next = double(load_data(myrecording.channels(1), nimg+1));

      [x,y,u,v] = matpiv_nfft(img, img_next, windows, 1/32, threshs, mask);
      empties = (u == 0 & v == 0);
      u(empties) = NaN;
      v(empties) = NaN;

      subplot(2,2,1);imagesc(img)
      subplot(2,2,2);imagesc(img_next)
      subplot(2,2,3);imagesc(mask);
      subplot(2,2,4);quiver(x,-y,u,-v)

      keyboard

      prev_indx = nimg;
    end
  end

  keyboard

  return;
end
