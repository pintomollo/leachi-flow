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
  %diff_thresh = 1;
  amp_thresh = 20;
  min_length = 1;
  min_branch = 10;
  prop_thresh = 0.3;
  nharmonics = 25;
  corr_thresh = 0.4;
  avg_thresh = 1e-5;
  pthresh = 1e-4;

  [nframes, img_size] = size_data(myrecording.channels(1));
  inelems = 1/(prod(img_size));

  if (isempty(segmentations.detections) || length(segmentations.detections) < nframes)
    detections = get_struct('detection', [1 nframes]);
  else
    detections = segmentations.detections(1:nframes);
  end

  vessel_width = ceil(max(b_leachi.vessel_width.mu) / opts.pixel_size);
  ampulla_width = b_leachi.ampulla.mu(1) / opts.pixel_size;
  %proj_dist = vessel_width * 0.75;
  proj_dist = 0.75;
  %sigma = min(b_leachi.blood_cell.mu(:,1)/(2*opts.pixel_size));
  sigma = min(b_leachi.blood_cell.mu(:,1)/(2*opts.pixel_size));

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
    % M-B-flow_1_1 ?
    %diff_thresh = 1 + 8./(1+exp(-((range(orig_img(:))/noise(2))-140)/12))
    % M-B-flow_1_10003 : 5.14
    diff_thresh = 1 + 8./(1+exp(-((range(orig_img(:))/noise(2))-130)/12));

    noise(1) = thresh;

    prev_img = gaussian_mex(orig_img, sigma);
    tmp_img = padarray(prev_img, [3 3]*vessel_width, NaN);
    %prev_mask = imdilate(imopen(tmp_img < noise(1) - diff_thresh*noise(2), disk2), disk2);
    %prev_mask = imdilate(imopen(imfill(tmp_img < noise(1) - amp_thresh*noise(2), 'holes'), disk2), disk2);
    prev_mask = imdilate(imopen(tmp_img < noise(1) - amp_thresh*noise(2), disk2), disk2);

    detections(1).noise = noise;

    if (opts.verbosity > 2)
      hfig=figure;
    end

    mask = zeros(img_size+6*vessel_width);
    for nimg=2:nframes
      new_img = double(load_data(myrecording.channels(1), nimg));
      new_img = gaussian_mex(new_img, sigma);

      [img_diff, moire] = immoire(new_img - prev_img, 5, 2.5*sigma);

      %orig_img = new_img;
      img = padarray(new_img, [3 3]*vessel_width, NaN);
      %curr_mask = imdilate(imopen(img < noise(1) - diff_thresh*noise(2), disk2), disk2);
      %curr_mask = imdilate(imopen(imfill(img < noise(1) - amp_thresh*noise(2), 'holes'), disk2), disk2);
      curr_mask = imdilate(imopen(img < noise(1) - amp_thresh*noise(2), disk2), disk2);

      img_diff = abs(padarray(img_diff, [3 3]*vessel_width, NaN));
      img_diff(prev_mask | curr_mask) = false;

      bw = img_diff > diff_thresh * noise(2);
      bw = bwareaopen(bw, ceil(5 / opts.pixel_size).^2);

      if (opts.verbosity > 2)
        %subplot(2,2,1);imagesc(prev_img)
        subplot(2,2,1);imagesc(prev_mask | curr_mask)
        subplot(2,2,2);imagesc(img_diff)
      end
      if (any(bw(:)))

        closed = imdilate(bw, disk1);
        open = imerode(closed, disk2);
        mask = mask + open;

        if (opts.verbosity > 2)
          subplot(2,2,3);imagesc(bw);
          subplot(2,2,4);imagesc(mask)
        end
        %props = sum(open(:)) * inelems;
        %title(props)
      end
      %title(nimg)
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
      if (opts.verbosity > 2)
        print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_masking.png'], '-dpng');
        delete(hfig)
      end
      close(hwait);
    end

    mask = imnorm(mask);
    %orig_mask = mask;

    %figure;subplot(2,2,1)
    %imagesc(mask)

    mask = (mask > prop_thresh);

    %subplot(2,2,2)
    %imagesc(mask)

    mask = bwareaopen(mask, vessel_width^2);

    %subplot(2,2,3)
    %imagesc(mask)

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

    %subplot(2,2,4)
    %imagesc(mask);hold on
    %plot(branches(:,1), branches(:,2), 'k');

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

    %subplot(2,2,4)
    %imagesc(mask);hold on
    %plot(branches(:,1), branches(:,2), 'k');

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
    %widths = 1/(proj_dist)^2;
    widths = 1 ./ (widths*proj_dist).^2;

    origin = [x1 y1].';
    params = [vects cross lens sqrt(lens .* widths)].';
  end

  tmp_img = double(load_data(myrecording.channels(1), 1));

  if (opts.verbosity > 2)
    hfig=figure;
    subplot(1,2,1);imagesc(tmp_img);hold on;plot(branches(:,1), branches(:,2), 'k');
    subplot(1,2,2);imagesc(mapping);hold on;plot(branches(:,1), branches(:,2), 'k');
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_mappings.png'], '-dpng');
    delete(hfig)
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
  snr = cell(nframes-1, 1);

  if (isempty(detections(2).carth) || all(isnan(detections(2).carth(:))))

    %figure;
    %%%%%%%%%%%%%%%%%%% SHOULD WORK ON THE DIFFERENCE BETWEEN FRAMES
    prev_img = double(load_data(myrecording.channels(1), 1));
    prev_img = gaussian_mex(prev_img, sigma);

    for nimg=2:nframes-1
      if (prev_indx == nimg-1)
        img = img_next;
        prev_diff = img_diff;
      else
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

      %[x,y,u,v,s] = matpiv_nfft(img, img_next, windows, 1/32, threshs, mask, 1.5);
      [x,y,u,v,s] = vessel_piv(prev_diff, img_diff, windows, 1/32, mapping, branches, threshs);
      %[x,y,u,v,s] = matpiv_nfft(prev_diff, img_diff, windows, 1/32, threshs, mask, 1);

      %for i=1:10
      %[x,y,u,v,s] = matpiv_nfft(guassian_mex(img, 0.67), gaussian_mex(img_next, 0.67), windows, 1/32, threshs, mask, i);

      %hold off;
      %imagesc(img_diff);
      %hold on;
      %quiver(x,y,u,v, 0);
      %title(nimg)
      %drawnow
      %end
      %keyboard

      %empties = (u == 0 & v == 0);
      %u(empties) = NaN;
      %v(empties) = NaN;
      %s(empties) = NaN;

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
      %subplot(2,2,4);imagesc((img_next-img) .* mask);hold on;
      %quiver(x,y,u,v, 'r');
      %hold off
      %drawnow
      %keyboard

      speed = bsxfun(@times, u(inside), params(1,:) .* sqrt(params(4,:))) + ...
               bsxfun(@times, v(inside), params(2,:) .* sqrt(params(4,:)));

      speed(others) = NaN;

      data{nimg} = speed*opts.pixel_size/opts.time_interval;
      detections(nimg).carth = [speed s(inside)];

      snr{nimg} = s(inside);

      prev_indx = nimg;
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
    %keyboard
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

      hfig=figure;
      subplot(1,3,1);plot(avgs);
      subplot(1,3,2);plot(avgs(:,new_empties));
      subplot(1,3,3);plot(avgs(:,~new_empties));
      print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_groupings.png'], '-dpng');
      delete(hfig);

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
    hfig = figure;hold on;
  end

  for i=1:size(avgs,2)
    if (opts.verbosity > 1)
      subplot(1,size(avgs,2), i); hold on;
      plot(pos*opts.time_interval, avgs(:,i), 'k');
    end

    avgs(:,i) = smooth(avgs(:,i), 0.2, 'rloess');
    %plot(smooth(avgs(:,i),15,'moving'), 'r');
    %plot(smooth(avgs(:,i),15,'sgolay'), 'g');
    %plot(smooth(avgs(:,i),0.2,'rlowess'), 'k');
    if (opts.verbosity > 1)
      plot(pos*opts.time_interval, avgs(:,i), 'Color', colors(i,:), 'LineWidth', 2);
    end
  end

  if (opts.verbosity > 1)
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
  %bads = cellfun('isempty', data);
  group_indxs = group_indxs*opts.time_interval;
  [gpos, indxi, indxj] = unique(group_indxs);

  %keyboard

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

  %{
  nbins = floor(max(SnR));
  colors = redbluemap(nbins+1);
  figure;hold on;
  for i=0:nbins
    goods = (SnR>i);
    scatter(group_indxs(goods), speeds(goods), 'MarkerEdgeColor', colors(i+1,:));
  end
  %}

  if (opts.verbosity > 1)
    hfig = figure;hold on;
    boxplot(speeds, group_indxs, 'position', gpos);
    set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');

    navgs = size(avgs, 2);
    colors = redbluemap(2*navgs);
    colors = colors(ceil(navgs/2)-1+[1:navgs],:);
    for i=1:navgs
      plot(gpos, avgs(:,i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_avg.png'], '-dpng');
    delete(hfig);
  end

  [t, v] = dp_flow(group_indxs, speeds, [500 128]);
  [pp] = csaps(t, v, 1/(t(end)));

  %{
  hfig = figure;hold on;
  scatter(group_indxs, speeds, 'k');
  fnplt(pp, 'r');
  print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_fits.png'], '-dpng');
  delete(hfig);
  %%
  %avg_val = mean(avgs, 2);

  nmax = 20;
  curr_speeds = speeds;
  colors=redbluemap(nmax);
  [mspeed, sspeed] = mymean(curr_speeds, 1, group_indxs);

  nbounds = max(range(gpos)*0.01, diff(gpos(1:2)));
  sval = smooth(sspeed, 0.05, 'rloess');
  thresh = median(sspeed);

  figure;hold on
  iters=[1:5:51];
  colors=redbluemap(length(iters));
  for n=1:length(iters)
    i = iters(n);
    t = imfilter(sspeed, ones(i,1)/i);
    plot(t, 'Color', colors(n,:))
  end
  plot([1 length(t)], [thresh thresh], 'k')

  keyboard

  w = (exp(-sspeed.^2 ./ (8*thresh^2)));
  %sval2 = sval .* (1-exp(-smooth(sspeed).^2 ./ (0.5*thresh^2)));
  sval = (sval .* (1-exp(-smooth(sspeed).^2 ./ (0.5*thresh^2)))) .* w + (1-w) .* sspeed;

  ngauss = length(sval)/25;
  h = fspecial('gaussian', [2*ceil(ngauss)+1 1], ngauss/3);
  gauss = imfilter(double(sval < thresh), h);

  %sval = sval .* (1-exp(-sspeed.^2 ./ (0.5*thresh^2)));
  bounds = (group_indxs <= gpos(1)+nbounds | group_indxs >= gpos(end)-nbounds);
  %all_goods = (sval(indxj) < thresh | bounds);
  all_goods = (gauss(indxj) > 0.5 | bounds);

  %all_goods2 = (sval2(indxj) < thresh | bounds);
  %all_goods3 = (sval3(indxj) < thresh | bounds);

  %[junk, pthresh] = csaps(group_indxs(all_goods), curr_speeds(all_goods));
  %[pp] = csaps(group_indxs(all_goods), curr_speeds(all_goods), pthresh/((opts.time_interval)^3));
  [pp] = csaps(group_indxs(all_goods), curr_speeds(all_goods), pthresh);
  avg_val = fnval(pp, gpos);

  hfig = figure;hold on;
  scatter(group_indxs, speeds, 'k');
  scatter(group_indxs(all_goods), speeds(all_goods), 'r');
  plot(gpos, avg_val, 'y');

  keyboard

  for i=1:nmax
    avg_speeds = avg_val(indxj);
    speed_dist = (curr_speeds - avg_speeds);
    %curr_thresh = median(speed_dist);

    %weight = 1-exp(-speed_dist.^2 / (2*(thresh^2)));
    %wspeed = (weight.*sign(speed_dist).*2*thresh) + avg_speeds;
    %[mspeed, sspeed] = mymean(curr_speeds, 1, group_indxs);
    %speed_dist = (curr_speeds - avg_speeds);
    %thresh = std(speed_dist);

    %goods = (sspeed < thresh & abs(mspeed - sign_val) < thresh);
    %all_goods = (goods(indxj));

    all_goods = ((abs(speed_dist) < thresh) & (sspeed(indxj) < thresh)) | bounds;

    %[junk, pthresh] = csaps(group_indxs(all_goods), curr_speeds(all_goods));
    %[pp] = csaps(group_indxs(all_goods), curr_speeds(all_goods), pthresh/((opts.time_interval)^3));
    [pp] = csaps(group_indxs(all_goods), curr_speeds(all_goods), pthresh);

    new_avg_val = fnval(pp, gpos);
    plot(gpos, new_avg_val, 'Color', colors(i,:));
    %plot(gpos, new_avg_val, 'r');

    dx_avg = (2*mean(abs(new_avg_val - avg_val)) / (range(avg_val) + range(new_avg_val)));
    avg_val = new_avg_val;

    if (dx_avg < avg_thresh)
      break;
    end
  end

  %if (opts.verbosity > 1)
    print(['-f' num2str(hfig)], ['./export/' myrecording.experiment '_fits.png'], '-dpng');
    delete(hfig);
  %end
  %}

  %{
  sign_val = mean(avgs, 2);
  %tmp_speeds = speeds - sign_val(indxj);
  %thresh = std(tmp_speeds)/2;
  %goods = (speeds < sign_val(indxj) + thresh & speeds > sign_val(indxj) - thresh) & ...
  %        (~isnan(group_indxs) & ~isnan(speeds));

  avg_speeds = sign_val(indxj);
  speed_dist = (speeds - avg_speeds);
  [mspeed, sspeed] = mymean(speeds, 1, group_indxs);

  thresh = std(speed_dist);
  weight = 1-exp(-speed_dist.^2 / (2*(thresh^2)));
  tmp_speed = (weight.*sign(speed_dist).*2*thresh) + avg_speeds;

  goods = (sspeed < thresh & abs(mspeed - sign_val) < thresh);
  all_goods = (goods(indxj));

  nmax = 20;
  thresh = thresh / nmax;


  %if (opts.verbosity > 2)
    figure;hold on;
    scatter(group_indxs, speeds, 'k');

    niter=5;
    colors=redbluemap(niter+1);
  %%% NEED TO ITERATE OVER A FEW VALUES (1/dx)
  for i=0:niter
    if i==0
      [pp, pthresh] = csaps(group_indxs(all_goods), tmp_speed(all_goods));
    else
      [pp] = csaps(group_indxs(all_goods), tmp_speed(all_goods), pthresh/((opts.time_interval)^i));
    end
    ys = fnval(pp, gpos);
    plot(gpos, ys, 'Color', colors(i+1,:));
  end
  %end
  %}

  %{
  keyboard

  prev_params = -Inf;
  for i=1:nmax
    %[period, ampls, phases] = lsqmultiharmonic(group_indxs(goods), speeds(goods));
    [period, ampls, phases] = lsqmultiharmonic(group_indxs, tmp_speed, nharmonics);
    nharm = length(ampls);

    %sign_val = zeros(size(gpos));
    %for j=1:nharm
    %  sign_val = sign_val + ampls(j)*cos((j*gpos/period)*2*pi + phases(j));
    %end
    %thresh = max(ampls)/2;
    sign_val = cossum(gpos, [ampls.'; (period.')./[1:nharm]; phases.']);

    if (opts.verbosity > 2)
      plot(gpos, sign_val, 'k');
    end

    avg_speeds = sign_val(indxj);

    speed_dist = (speeds - avg_speeds);
    thresh = std(speed_dist);
    weight = 1-exp(-speed_dist.^2 / (2*(thresh^2)));
    tmp_speed = (weight.*sign(speed_dist).*2*thresh) + avg_speeds;

    %weight = exp(-(speeds - avg_speeds).^2 / (2*(((nmax-i)*thresh)^2)));
    %tmp_speed = avg_speeds.*(1-weight) + speeds.*weight;

    %goods = (speeds < sign_val(indxj) + thresh & speeds > sign_val(indxj) - thresh);

    dx = sum(abs([period ampls(1) phases(1)] - prev_params));
    prev_params = [period ampls(1) phases(1)];

    if (dx < 1e-3)
      break;
    end
  end
  %}

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
  %res.cluster = weight;
  %res.properties = [ampls(:).'; period*(1./[1:length(ampls)]); phases(:).'];
  res.cluster = [];
  res.properties = pp;

  if (opts.verbosity > 1)
    %[gpos2] = unique(group_indxs(goods));
    hfig = figure;hold on;
    boxplot(speeds, group_indxs, 'position', gpos);
    set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');

    plot(gpos, avg_val, 'k', 'LineWidth', 2);
    tmp = range(avg_val);
    ylim([-tmp tmp]);
    %tmp = 2*sum(ampls(:));
    %ylim([-tmp tmp]);
    %title(num2str(res.properties));
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
