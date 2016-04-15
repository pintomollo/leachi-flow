function [myrecording, opts] = leachi_ampullae(myrecording, opts)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

  if (isnumeric(myrecording))
    all_props = myrecording;
    b_leachi = get_struct('botrylloides_leachi');
    ampulla_width = b_leachi.ampulla.mu(1) / opts.pixel_size;
    frac_thresh = 0.95;
    nframes = max(all_props(:,end));
  else

  if (~isstruct(myrecording))
    if (ischar(myrecording))

      if (any(myrecording=='*'))
        files = dir(myrecording);
        [file_path, filename, ext] = fileparts(myrecording);
        for i=1:length(files)
          try
            leachi_ampullae(fullfile(file_path, files(i).name), opts, true);
          catch
            disp(lasterr)
          end
        end

        myrecording = [];
        opts = [];

        return;
      end

      [file_path, filename, ext] = fileparts(myrecording);
      if (~strncmp(ext, '.tif', 4) && ~strncmp(ext, '.tiff', 5))
        myrecording = convert_movie(myrecording, false);
      end

      [myrecording, opts] = inspect_recording(myrecording, opts, batch_mode);
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

  [nframes, img_size] = size_data(myrecording.channels(1));

  if (isempty(segmentations.detections) || length(segmentations.detections) < nframes)
    detections = get_struct('detection', [1 nframes]);
  else
    detections = segmentations.detections(1:nframes);
  end

  if (~isempty(detections(1).cluster))
    mapping = detections(1).cluster;
    noise = detections(1).noise;
    mask = ~logical(mapping);
  else
    mask = true(img_size);
    noise = [];
  end

  b_leachi = get_struct('botrylloides_leachi');
  ampulla_width = b_leachi.ampulla.mu(1) / opts.pixel_size;
  amp_thresh = 20;
  disk1 = strel('disk', ceil(ampulla_width/2));
  disk2 = strel('disk', ceil(ampulla_width/4));
  nmagic = ceil(nframes/2);
  frac_thresh = 0.95;

  if (opts.verbosity > 1)
    hwait = waitbar(0,'Detecting ampullae...','Name','B. leachi ampullae');
  end

  if (isempty(noise))
    orig_img = double(load_data(myrecording.channels(1), 1));

    bkg = gaussian_mex(orig_img, ampulla_width);
    thresh = opthr(bkg);

    noise = estimate_noise(orig_img);
    diff_thresh = 1 + 8./(1+exp(-((range(orig_img(:))/noise(2))-130)/12));

    noise(1) = thresh;
  end

  if (opts.verbosity > 1)
    hfig=figure;
    hfig = hfig.Number;
    colormap(hfig, redbluemap);
  end

  ampullae = cell(nframes, 1);
  nampullae = zeros(nframes, 1);
  for nimg=1:nframes
    img = double(load_data(myrecording.channels(1), nimg));
    ampulla = imdilate(imopen(img < noise(1) - amp_thresh*noise(2), disk2), disk2);
    ampulla_props = regionprops(ampulla, 'Centroid', 'Area');

    nampullae(nimg) = length(ampulla_props);
    props = NaN(nampullae(nimg), 4);
    for j=1:nampullae(nimg)
      props(j, 1) = ampulla_props(j).Area;
      props(j, 2:3) = ampulla_props(j).Centroid(:).';
      props(j, 4) = nimg;
    end

    ampullae{nimg} = props;
    %detections(nimg).cluster = ampullae;

    if (opts.verbosity > 2 || (opts.verbosity>1 && nimg==nmagic))
      h=subplot(2,1,1,'Parent',hfig);imagesc(img,'Parent',h);
      h=subplot(2,1,2,'Parent',hfig);imagesc(ampulla,'Parent',h)
    end

    if (opts.verbosity > 1)
      waitbar(nimg/nframes,hwait);
    end
    if (opts.verbosity > 2)
      drawnow
    end
  end

  if (opts.verbosity > 1)
    close(hwait);
  end

  all_props = cat(1, ampullae{:});
  end

  clusts = simple_clustering(all_props, ampulla_width);
  ids = unique(clusts);

  nclusts = length(ids);
  clusts_lengths = zeros(nclusts, 1);

  for i=1:nclusts
    clusts_lengths(i) = length(unique(all_props(clusts==i, end)));
  end

  goods = (clusts_lengths / nframes) > frac_thresh;
  ids = ids(goods);

  goods = ismember(clusts, ids);
  clusts = clusts(goods);
  all_props = all_props(goods, :);
  nclusts = length(ids);

  all_ampullae = NaN(nclusts, nframes);
  for i=1:nclusts
    tmp_pts = all_props(clusts==ids(i),:);
    [frames, indxi, indxj] = unique(tmp_pts(:,end));

    all_ampullae(i, frames) = tmp_pts(indxi, 1);

    if (length(indxi) ~= length(indxj))
      tmp_pts(indxi,:) = [];

      for j=1:size(tmp_pts, 1)
        [frames, indxi, indxj] = unique(tmp_pts(:,end));

        all_ampullae(i, frames) = all_ampullae(i, frames) + tmp_pts(indxi, 1);
        tmp_pts(indxi,:) = [];

        if (isempty(tmp_pts))
          break;
        end
      end
    end
  end

  %%smooth

  keyboard

  return;
end

function clusters = simple_clustering(pts, thresh)

  npts = size(pts, 1);
  clusters = zeros(npts, 1);
  count = 1;
  thresh = thresh^2;

  for i=1:npts
    remainings = (~clusters);

    if (~any(remainings))
      break;
    end

    if (remainings(i))
      tmp_pts = pts(remainings, 2:3);

      new = logical(clusters(remainings));
      new(1) = true;

      for j=1:length(new)
        dist = bsxfun(@minus, tmp_pts(~new,1), tmp_pts(new,1).').^2 + bsxfun(@minus, tmp_pts(~new,2), tmp_pts(new,2).').^2;
        goods = any(dist <= thresh, 2);

        if (~any(goods))
          break;
        end

        new(~new) = goods;
      end

      new = new * count;
      clusters(remainings) = new;
      count = count + 1;
    end
  end

  return;
end
