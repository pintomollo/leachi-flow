function [myrecording, opts] = leachi_ampullae(myrecording, opts)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

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

  noise = detections(1).noise;

  b_leachi = get_struct('botrylloides_leachi');
  ampulla_width = b_leachi.ampulla.mu(1) / opts.pixel_size;
  amp_thresh = 20;
  disk1 = strel('disk', ceil(ampulla_width/2));
  disk2 = strel('disk', ceil(ampulla_width/4));
  nmagic = ceil(nframes/2);
  frac_thresh = 0.95;
  dt = opts.time_interval;

  if (isempty(noise))
    orig_img = double(load_data(myrecording.channels(1), 1));

    bkg = gaussian_mex(orig_img, ampulla_width);
    thresh = opthr(bkg);

    noise = estimate_noise(orig_img);
    diff_thresh = 1 + 8./(1+exp(-((range(orig_img(:))/noise(2))-130)/12));

    noise(1) = thresh;
    detections(1).noise = noise;
  end

  if (isempty(detections(1).cluster) || any(isnan(detections(1).cluster(:))))

    if (opts.verbosity > 1)
      hwait = waitbar(0,'Detecting ampullae...','Name','B. leachi ampullae');
      hfig=figure;
      hfig = hfig.Number;
      colormap(hfig, redbluemap);
    end

    ampullae = cell(nframes, 1);
    nampullae = zeros(nframes, 1);
    for nimg=1:nframes
      img = double(load_data(myrecording.channels(1), nimg));
      %ampulla = imdilate(imopen(img < noise(1) - amp_thresh*noise(2), disk2), disk2);
      ampulla = imclose(imopen(img < noise(1) - amp_thresh*noise(2), disk2), disk2);
      ampulla_props = regionprops(ampulla, 'Centroid', 'Area');

      nampullae(nimg) = length(ampulla_props);
      props = NaN(nampullae(nimg), 4);
      for j=1:nampullae(nimg)
        props(j, 1) = ampulla_props(j).Area;
        props(j, 2:3) = ampulla_props(j).Centroid(:).';
        props(j, 4) = nimg;
      end

      ampullae{nimg} = props;
      detections(nimg).cluster = ampullae{nimg};

      if (opts.verbosity > 2 || (opts.verbosity>1 && nimg==nmagic))
        h=subplot(2,2,1,'Parent',hfig);imagesc(img,'Parent',h);
        h=subplot(2,2,2,'Parent',hfig);imagesc(ampulla,'Parent',h)
        h=subplot(2,2,3,'Parent',hfig);imagesc(bwlabel(ampulla),'Parent',h)
        h=subplot(2,2,3,'Parent',hfig);set(h, 'XLim', [1 img_size(2)], 'YLim', [1 img_size(1)], 'YDir', 'reverse', 'NextPlot', 'add');
        scatter(h, props(:,2), props(:,3), props(:,1)*get_pts_resolution(h)^2);
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

    segmentations.detections = detections;
    myrecording.segmentations = segmentations;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  else
    ampullae = cell(nframes, 1);
    for nimg=1:nframes
      ampullae{nimg} = detections(nimg).cluster;
    end
  end

  all_props = cat(1, ampullae{:});

  clusts = simple_clustering(all_props, ampulla_width);
  avg_sizes = mymean(all_props(:,1), 1, clusts);
  ids = unique(clusts);

  nclusts = length(ids);
  clusts_lengths = zeros(nclusts, 1);

  for i=1:nclusts
    clusts_lengths(i) = length(unique(all_props(clusts==i, end)));
  end

  goods = (clusts_lengths / nframes) > frac_thresh;
  ids = ids(goods);
  avg_sizes = avg_sizes(goods);
  [junk, col_ind] = sort(avg_sizes);
  [junk, col_ind] = sort(col_ind);

  goods = ismember(clusts, ids);
  clusts = clusts(goods);
  contam = all_props(~goods, :);
  all_props = all_props(goods, :);
  nclusts = length(ids);

  colors = parula(2*nclusts);

  tmp_indxs = ones(max(ids), 1);
  tmp_indxs(ids) = col_ind;
  all_colors = colors(2*tmp_indxs(clusts),:);

  figure;h=axes();
  set(h, 'XLim', [1 img_size(2)], 'YLim', [1 img_size(1)], 'YDir', 'reverse', 'NextPlot', 'add');
  scatter(all_props(1:10:end,2), all_props(1:10:end,3), all_props(1:10:end,1)*get_pts_resolution(h)^2, all_colors(1:10:end,:));
  scatter(contam(:,2), contam(:,3), contam(:,1)*get_pts_resolution(h)^2, 'k');

  all_ampullae = NaN(nclusts, nframes);
  %all_fits = cell(nclusts, 1);

  figure;hold on;

  all_fits = NaN(nclusts, nframes);
  all_extrema = cell(nclusts, 2);
  t = [1:nframes]*dt;
  for i=1:nclusts
    tmp_pts = all_props(clusts==ids(i),:);
    [frames, indxi, indxj] = unique(tmp_pts(:,end));

    all_ampullae(i, frames) = tmp_pts(indxi, 1);

    if (length(indxi) ~= length(indxj))
      tmp_pts(indxi,:) = [];

      for j=1:size(tmp_pts, 1)
        [frames, indxi, indxj] = unique(tmp_pts(:,end));

        all_ampullae(i, frames) = all_ampullae(i, frames) + tmp_pts(indxi, 1).';
        tmp_pts(indxi,:) = [];

        if (isempty(tmp_pts))
          break;
        end
      end
    end

    goods = ~isnan(all_ampullae(i,:));
    [pp] = csaps(t(goods), all_ampullae(i,goods), 1/(10*t(end)));
    all_fits(i,:) = fnval(pp, t);
    [xmax,imax,xmin,imin] = local_extrema(all_fits(i,:));
    all_extrema{i,1} = [xmax, imax];
    all_extrema{i,2} = [xmin, imin];

    plot(t, all_ampullae(i,:), 'Color', colors(2*col_ind(i),:));
    plot(t, all_fits(i,:), 'Color', colors(2*col_ind(i)-1,:));
    scatter(t(all_extrema{i,1}(:,2)), 500+all_extrema{i,1}(:,1), 'v', 'MarkerEdgeColor', colors(2*col_ind(i)-1,:))
    scatter(t(all_extrema{i,2}(:,2)), -500+all_extrema{i,2}(:,1), '^', 'MarkerEdgeColor', colors(2*col_ind(i)-1,:))
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
