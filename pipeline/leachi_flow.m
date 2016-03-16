function [myrecording, opts] = leachi_flow(myrecording, opts, batch_mode)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

  if (nargin < 3)
    batch_mode = false;
  end

  if (~isstruct(myrecording))
    if (ischar(myrecording))

      if (any(myrecording=='*'))
        files = dir(myrecording);
        [file_path, filename, ext] = fileparts(myrecording);
        for i=1:length(files)
          try
            leachi_flow(fullfile(file_path, files(i).name), opts, true);
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

  b_leachi = get_struct('botrylloides_leachi');

  diff_thresh = 6;
  amp_thresh = 20;
  min_length = 1;
  min_branch = 10;
  prop_thresh = 0.3;
  nharmonics = 25;
  corr_thresh = 0.4;
  avg_thresh = 1e-5;
  pthresh = 1e-4;

  vessel_width = ceil(max(b_leachi.vessel_width.mu) / opts.pixel_size);
  ampulla_width = b_leachi.ampulla.mu(1) / opts.pixel_size;
  proj_dist = 0.75;
  sigma = min(b_leachi.blood_cell.mu(:,1)/(2*opts.pixel_size));

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));
  nmagic = ceil(nframes/2);
  magic_img = [];

  if (isempty(segmentations.detections) || length(segmentations.detections) < nframes)
    detections = get_struct('detection', [1 nframes]);
  else
    detections = segmentations.detections(1:nframes);
  end

  if (isempty(detections(1).cluster))
    if (opts.verbosity > 1)
      hwait = waitbar(0,'Computing mask...','Name','B. leachi flow');
    end

    orig_img = double(load_data(myrecording.channels(1), 1));

    disk1 = strel('disk', 2*vessel_width);
    disk2 = strel('disk', vessel_width);

    bkg = gaussian_mex(orig_img, ampulla_width);
    thresh = opthr(bkg);

    noise = estimate_noise(orig_img);
    diff_thresh = 1 + 8./(1+exp(-((range(orig_img(:))/noise(2))-130)/12));

    noise(1) = thresh;

    prev_img = gaussian_mex(orig_img, sigma);
    tmp_img = padarray(prev_img, [3 3]*vessel_width, NaN);
    prev_mask = imdilate(imopen(tmp_img < noise(1) - amp_thresh*noise(2), disk2), disk2);

    detections(1).noise = noise;

    if (opts.verbosity > 1)
      hfig=figure;
      colormap(hfig, redbluemap);
    end

    mask = zeros(img_size+6*vessel_width);
    for nimg=2:nframes
      new_img = double(load_data(myrecording.channels(1), nimg));
      new_img = gaussian_mex(new_img, sigma);

      [img_diff, moire] = immoire(new_img - prev_img, 5, 2.5*sigma);
      if (nimg==nmagic)
        magic_img = new_img;
      end

      img = padarray(new_img, [3 3]*vessel_width, NaN);
      curr_mask = imdilate(imopen(img < noise(1) - amp_thresh*noise(2), disk2), disk2);

      img_diff = abs(padarray(img_diff, [3 3]*vessel_width, NaN));
      img_diff(prev_mask | curr_mask) = false;

      bw = img_diff > diff_thresh * noise(2);
      bw = bwareaopen(bw, ceil(5 / opts.pixel_size).^2);

      if (opts.verbosity > 2 || (opts.verbosity>1 && nimg==nmagic))
        h=subplot(2,2,1,'Parent',hfig);imagesc(magic_img,'Parent',h);
        h=subplot(2,2,2,'Parent',hfig);imagesc(img_diff,'Parent',h)
      end
      if (any(bw(:)))

        closed = imdilate(bw, disk1);
        open = imerode(closed, disk2);
        mask = mask + open;

        if (opts.verbosity > 2 || (opts.verbosity>1 && nimg==nmagic))
          h=subplot(2,2,3,'Parent',hfig);imagesc(bw,'Parent',h);
          h=subplot(2,2,4,'Parent',hfig);imagesc(mask,'Parent',h)
        end
      end
      detections(nimg).carth = [];

      prev_img = new_img;
      prev_mask = curr_mask;

      if (opts.verbosity > 1)
        waitbar(nimg/nframes,hwait);
      end
      if (opts.verbosity > 2)
        drawnow
      end
    end

    if (opts.verbosity > 1)
      plot2svg(['./export/SVG/' myrecording.experiment '_maskinga.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_maskinga.png'], '-dpng');

      colormap(hfig, gray);
      plot2svg(['./export/SVG/' myrecording.experiment '_maskingb.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_maskingb.png'], '-dpng');
      delete(hfig)

      close(hwait);
    end

    mask = imnorm(mask);

    if (opts.verbosity > 1)
      hfig=figure;
      colormap(hfig, redbluemap);
      h=subplot(2,2,1,'Parent',hfig);imagesc(mask,'Parent',h)
    end

    mask = (mask > prop_thresh);

    if (opts.verbosity > 1)
      h=subplot(2,2,2,'Parent',hfig);imagesc(mask,'Parent',h)
    end

    mask = bwareaopen(mask, vessel_width^2);

    if (opts.verbosity > 1)
      h=subplot(2,2,3,'Parent',hfig);imagesc(mask,'Parent',h)
    end

    dist = bwdist(~mask);
    mask = bwmorph(mask, 'thin', Inf);
    dist = dist(mask);

    mask = mask(3*vessel_width+[1:img_size(1)], 3*vessel_width+[1:img_size(2)]);
    [icoord, jcoord] = find(mask);
    centers = [jcoord, icoord];

    branches = sort_shape([centers dist], min_branch, vessel_width / 2);
    if (isempty(branches))
      error('nothing');
    end
    branches(:,4) = NaN;

    detections(1).properties = branches;

    if (opts.verbosity > 1)
      h=subplot(2,2,4,'Parent',hfig);imagesc(magic_img,'Parent',h);
      show_vessels(branches, h);

      plot2svg(['./export/SVG/' myrecording.experiment '_branchingsa.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_branchingsa.png'], '-dpng');

      colormap(hfig, gray)
      plot2svg(['./export/SVG/' myrecording.experiment '_branchingsb.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_branchingsb.png'], '-dpng');
      delete(hfig)
    end

    x1 = branches(1:3:end, 1);
    x2 = branches(2:3:end, 1);
    y1 = branches(1:3:end, 2);
    y2 = branches(2:3:end, 2);
    widths = branches(1:3:end, 3);
    widths(widths<vessel_width) = vessel_width;
    all_signs = branches(1:3:end, 4);

    branches_size = length(x1);

    vects = [x2-x1 y2-y1];
    lens = 1 ./ sum(vects.^2, 2);
    cross = (x2.*y1 - y2.*x1);
    widths = 1./((widths*proj_dist).^2);

    origin = [x1 y1].';
    params = [vects cross lens sqrt(lens .* widths)].';

    [X,Y] = meshgrid([1:img_size(2)], [1:img_size(1)]);

    dists = bsxfun(@times, (bsxfun(@plus, ...
                              bsxfun(@times, X(:), params(2, :)) - ...
                                bsxfun(@times, Y(:), params(1,:)), ...
                              params(3,:))).^2, ...
                           params(4,:) .* widths.');

    nodes = bsxfun(@times, (bsxfun(@minus, X(:), [x1;x2].').^2 + ...
                            bsxfun(@minus, Y(:), [y1;y2].').^2), ...
                            [widths; widths].');

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

    x1 = branches(1:3:end, 1);
    x2 = branches(2:3:end, 1);
    y1 = branches(1:3:end, 2);
    y2 = branches(2:3:end, 2);
    widths = branches(1:3:end, 3);
    widths(widths<vessel_width) = vessel_width;
    branches_size = length(x1);

    if (size(branches, 2)>3)
      all_signs = branches(1:3:end, 4);
    else
      all_signs = NaN(branches_size, 1);
    end

    vects = [x2-x1 y2-y1];
    lens = 1 ./ sum(vects.^2, 2);
    cross = (x2.*y1 - y2.*x1);
    widths = 1 ./ (widths*proj_dist).^2;

    origin = [x1 y1].';
    params = [vects cross lens sqrt(lens .* widths)].';
  end

  if (isempty(magic_img))
    magic_img = double(load_data(myrecording.channels(1), nmagic));
  end

  if (opts.verbosity > 1)
    hfig=figure;
    h=subplot(1,2,1,'Parent',hfig);imagesc(magic_img,'Parent',h);
    show_vessels(branches, h);
    h=subplot(1,2,2,'Parent',hfig);imagesc(mapping,'Parent',h);

    colormap(h, gray);
    plot2svg(['./export/SVG/' myrecording.experiment '_mappingsa.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_mappingsa.png'], '-dpng');

    colormap(h, brewermap([], 'Purples'));
    plot2svg(['./export/SVG/' myrecording.experiment '_mappingsb.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_mappingsb.png'], '-dpng');

    delete(hfig)
  end

  prev_indx = -1;
  real_mapping = [];

  windows = vessel_width * [4 3 2 1 1].';
  threshs = [5; 5; 3; 3];

  data = cell(nframes-1, 1);
  snr = cell(nframes-1, 1);

  if (isempty(detections(2).carth) || all(isnan(detections(2).carth(:))))

    if (opts.verbosity > 1)
      hfig=figure;
      colormap(hfig, redbluemap);
    end

    prev_img = double(load_data(myrecording.channels(1), 1));
    prev_img = gaussian_mex(prev_img, sigma);

    for nimg=2:nframes-1
      if (prev_indx == nimg-1)
        img = img_next;
        prev_diff = img_diff;
      else
        prev_img = double(load_data(myrecording.channels(1), nimg-1));
        img = double(load_data(myrecording.channels(1), nimg));
        img = gaussian_mex(img, sigma);
        [prev_diff, moire] = immoire(img - prev_img, 5, 2.5*sigma);
        prev_diff = imdetrend(prev_diff);
      end
      img_next = double(load_data(myrecording.channels(1), nimg+1));
      img_next = gaussian_mex(img_next, sigma);
      [img_diff, moire] = immoire(img_next - img, 5, 2.5*sigma);
      img_diff = imdetrend(img_diff);

      %%%%%%% COULD FILTER OUT VECTORS THAT ARE NOT // WITH THE CENTERS. EITHER DURING OR AFTER THE PIV

      [x,y,u,v,s] = vessel_piv(prev_diff, img_diff, windows, 1/32, mapping, branches, threshs);

      if (isempty(real_mapping))
        tmp_vals = bilinear_mex(mapping, x, y);
        real_mapping = blockproc(tmp_vals, [1 1], @local_mapping, 'BorderSize', [1 1], 'TrimBorder', false);

        inside = (real_mapping > 0);

        ngoods = sum(inside(:));
        subs = sub2ind([ngoods branches_size], [1:ngoods].', real_mapping(inside));

        others = ~ismember([1:ngoods*branches_size], subs);
        one = ones(1, branches_size);
      end

      if (opts.verbosity > 2 || (opts.verbosity>1 && nimg==nmagic))
        h=subplot(2,2,1,'Parent',hfig);imagesc(prev_diff,'Parent',h);
        h=subplot(2,2,2,'Parent',hfig);imagesc(img_diff,'Parent',h);
        h=subplot(2,2,3,'Parent',hfig);imagesc(mapping,'Parent',h);
        h=subplot(2,2,4,'Parent',hfig, 'NextPlot', 'add');imagesc(img,'Parent',h);
        quiver(h, x,y,u,v, 'r');

        if (opts.verbosity > 2)
          drawnow
        end
      end

      speed = bsxfun(@times, u(inside), params(1,:) .* sqrt(params(4,:))) + ...
               bsxfun(@times, v(inside), params(2,:) .* sqrt(params(4,:)));

      speed(others) = NaN;

      data{nimg} = speed*opts.pixel_size/opts.time_interval;
      detections(nimg).carth = [speed s(inside)];

      snr{nimg} = s(inside);

      prev_indx = nimg;
    end

    if (opts.verbosity > 1)
      plot2svg(['./export/SVG/' myrecording.experiment '_piva.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_maskinga.png'], '-dpng');

      colormap(hfig, gray);
      plot2svg(['./export/SVG/' myrecording.experiment '_pivb.svg'], hfig, 'png');
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_maskingb.png'], '-dpng');
      delete(hfig)
    end

    segmentations.detections = detections;
    myrecording.segmentations = segmentations;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  else
    for nimg=2:nframes-1
      data{nimg} = detections(nimg).carth(:,1:end-1)*opts.pixel_size/opts.time_interval;
      snr{nimg} = detections(nimg).carth(:,end);
    end
  end

  data = data(2:end);
  snr = snr(2:end);

  ndata = length(data);
  avgs = cellfun(@nanmedian, data, 'UniformOutput', false);
  avgs = cat(1, avgs{:});
  pos = [1:ndata];
  content = (sum(isnan(avgs), 1) < ndata/20);
  nils = any(isnan(avgs(:,content)), 2);
  pos = pos(~nils);
  avgs = avgs(~nils,:);

  if (all(isnan(all_signs)) || all(all_signs == 0))
    if (sum(content) == 1)
      all_signs = content;
      empty_branches = (all_signs == 0);
      sames = all_signs(~empty_branches);
      sames = sames(:).';
    else
      tmp_avgs = avgs;
      tmp_avgs(isnan(tmp_avgs)) = 0;
      corrs = corr(tmp_avgs);

      C = corrs - eye(branches_size);
      sC = (C >= 0);
      aC = abs(C);
      aC(aC < corr_thresh) = NaN;
      groups = NaN(branches_size);
      groups(:,1) = [1:branches_size];

      empty_branches = all(isnan(C), 1);

      for i=1:branches_size^2
        [val, indxi] = max(aC,[],1);
        [val, indxj] = max(val,[],2);

        if (val == 0)
          tmp_groups = find(~isnan(groups(:,2)));
          if (length(tmp_groups) > 1)
            aC(indxi,indxj) = abs(C(indxi,indxj));
            [val, indxi] = max(aC,[],1);
            [val, indxj] = max(val,[],2);
          else
            break;
          end
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

      new_empties = xor(isnan(groups(:,1)), isnan(groups(:,2)));
      values = groups(~isnan(groups));
      [junk, indx] = sort(abs(values));
      sames = sign(values(indx));

      if (opts.verbosity > 1)
        hfig=figure;
        h=subplot(1,3,1,'Parent',hfig);plot(avgs,'Parent',h);
        h=subplot(1,3,2,'Parent',hfig);plot(avgs(:,new_empties),'Parent',h);
        h=subplot(1,3,3,'Parent',hfig);plot(avgs(:,~new_empties),'Parent',h);
        plot2svg(['./export/SVG/' myrecording.experiment '_groupings.svg'], hfig, 'png');
        print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_groupings.png'], '-dpng');
        delete(hfig);
      end

      empty_branches(new_empties) = true;

      sames = sames(~empty_branches);
      sames = sames(:).';

      all_signs = double(~empty_branches);
      all_signs(~empty_branches) = sames;
    end

    branches(1:3:end,4) = all_signs(:);
    detections(1).properties = branches;

    segmentations.detections = detections;
    myrecording.segmentations = segmentations;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  else
    empty_branches = (all_signs == 0);
    sames = all_signs(~empty_branches);

    sames = sames(:).';
  end

  avgs = bsxfun(@times, avgs(:,~empty_branches), sames);

  if (opts.verbosity > 1)
    navgs = size(avgs, 2);
    colors = redbluemap(2*navgs);
    colors = colors(ceil(navgs/2)-1+[1:navgs],:);
    hfig = figure;
  end

  for i=1:size(avgs,2)
    if (opts.verbosity > 1)
      h=subplot(1,size(avgs,2), i);
      set(h, 'NextPlot', 'add');
      plot(pos*opts.time_interval, avgs(:,i), 'k', 'Parent', h);
    end

    avgs(:,i) = smooth(avgs(:,i), 0.2, 'rloess');
    if (opts.verbosity > 1)
      plot(pos*opts.time_interval, avgs(:,i), 'Color', colors(i,:), 'LineWidth', 2, 'Parent', h);
    end
  end

  if (opts.verbosity > 1)
    plot2svg(['./export/SVG/' myrecording.experiment '_branches.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_branches.png'], '-dpng');
    delete(hfig);
  end

  speeds = [];
  group_indxs = [];
  avgs_indxs = [];
  SnR = [];
  for i = pos
    tmp_all = bsxfun(@times, data{i}(:,~empty_branches), sames);
    s = snr{i};
    for j=1:size(avgs,2)
      tmp = tmp_all(:,j);
      SnR = [SnR; s(isfinite(tmp))];
      tmp = tmp(isfinite(tmp));
      speeds = [speeds; tmp];
      group_indxs = [group_indxs; ones(size(tmp))*i];
      avgs_indxs = [avgs_indxs; ones(size(tmp))*j];
    end
  end
  group_indxs = group_indxs*opts.time_interval;
  [gpos, indxi, indxj] = unique(group_indxs);

  if (opts.verbosity > 1)
    hfig = figure;
    haxes = axes('Parent', hfig, 'NextPlot', 'add');
    boxplot(haxes, speeds, group_indxs, 'position', gpos);
    set(haxes, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');

    navgs = size(avgs, 2);
    colors = redbluemap(2*navgs);
    colors = colors(ceil(navgs/2)-1+[1:navgs],:);
    for i=1:navgs
      plot(gpos, avgs(:,i), 'Color', colors(i,:), 'LineWidth', 2, 'Parent', haxes);
    end
    plot2svg(['./export/SVG/' myrecording.experiment '_avg.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_avg.png'], '-dpng');
    delete(hfig);
  end

  [t, v] = dp_flow(group_indxs, speeds, [0.35 6]);

  [pp] = csaps(t, v, 1/(t(end)));
  avg_val = fnval(pp, gpos);

  if (opts.verbosity > 1)
    hfig = figure;
    haxes = axes('Parent', hfig, 'NextPlot', 'add');
    scatter(haxes, group_indxs, speeds, 'k');
    plot(t, v, 'y', 'Parent', haxes)
    plot(gpos, avg_val, 'r', 'Parent', haxes);

    plot2svg(['./export/SVG/' myrecording.experiment '_fits.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_fits.png'], '-dpng');
    delete(hfig);
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
  res.cluster = [];
  res.properties = pp;

  if (opts.verbosity > 1)
    hfig = figure;
    haxes = axes('Parent', hfig, 'NextPlot', 'add');
    boxplot(haxes, speeds, group_indxs, 'position', gpos);
    set(haxes, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');
    plot(gpos, avg_val, 'k', 'LineWidth', 2, 'Parent', haxes);
    tmp = range(avg_val);
    ylim(haxes, [-tmp tmp]);

    plot2svg(['./export/SVG/' myrecording.experiment '_f.svg'], hfig, 'png');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_f.png'], '-dpng');
    delete(hfig);
  end

  trackings.detections = res;

  myrecording.segmentations = segmentations;
  myrecording.trackings = trackings;

  save([myrecording.experiment '.mat'], 'myrecording', 'opts');

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
