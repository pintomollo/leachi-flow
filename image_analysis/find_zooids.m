function [all_coords, sizes] = find_zooids(img, systems, params)

  sizes = [];
  if (nargin == 1 && isstruct(img))
    mycolony = img;
    nchannels = length(mycolony.channels);
    sizes = zeros(1, nchannels);

    fprintf('Locating zooids :     ');

    for i=1:nchannels
      fprintf('\b\b\b%3d', i);

      img = imread(mycolony.channels(i).fname);
      if (mycolony.channels(i).normalize)
        img = imnorm(img);
      end

      params = [8/mycolony.channels(i).pixel_size mycolony.channels(i).amplitude];
      systems = mycolony.channels(i).system;

      zooids = find_zooids(img, systems, params);
      sizes(i) = size(zooids, 1);

      mycolony.channels(i).zooids = zooids;
    end

    fprintf('\b\b\b\bdone\n');
    fprintf('Total found in %s : %d\n', mycolony.experiment, sum(sizes));

    all_coords = mycolony;

    return;
  end

  %orig_img = img;
  %orig_sys = systems;
  if (isempty(img) || isempty(systems))
    all_coords = NaN(1,2);

    return;
  end

  if (size(img, 3) > 1)
    img = rgb2gray(img);
  end

  rads = 3;
  ampl = -1;

  if (nargin > 2)
    rads = params(1);

    if (numel(params) > 1)
      ampl = params(2);
    end
  end

  pext = 5*rads;

  vects = [1 1; diff(systems)];
  dist = sqrt(sum(vects.^2, 2));
  vects = bsxfun(@rdivide, vects, dist);

  goods = (dist~=0);

  vects = vects(goods, :);
  dist = dist(goods);
  systems = systems(goods, :);

  gaps = all(isnan(systems), 2);

  first = [true; gaps(1:end-1)];
  second = [false; first(1:end-1)];

  systems(first,:) = systems(second,:) + bsxfun(@times, -vects(second,:), dist(second)+pext);

  last = [gaps(2:end); false];
  prev = [last(2:end); false];

  systems(last,:) = systems(prev,:) + bsxfun(@times, vects(last,:), dist(last)+pext);

  [perp, paths, dpos] = perpendicular_sampling(img, systems);

  if (~iscell(perp))
    perp = {perp};
    paths = {paths};
  end

  weight = exp(-(dpos.^2)/(2*(2*rads)^2));
  all_coords = NaN(0,5);

  for i=1:length(perp)

    proj = perp{i};
    pos = paths{i};
    pos = reshape(pos, [size(proj) 2]);

    if (rads>0)
      proj = gaussian_mex(double(proj), rads);
    else
      proj = double(proj);
    end

    evol = proj(:,(end-1)/2 + 1);
    [xmax, imax, xmin, imin] = local_extrema(evol, rads);

    center_proj = 1-bsxfun(@times, 1-imnorm(proj), weight);

    nzooids = length(imin);
    coords = NaN(nzooids, 5);

    prev_min = NaN;
    prev_pos = 0;
    for j=1:length(imin)
      indxi = imin(j);
      curr_row = center_proj(indxi,:);
      [val, indxj] = min(curr_row);

      coords(j,1:2) = squeeze(pos(indxi, indxj, :));

      min_val = proj(indxi, indxj);
      max_val = max([xmax(imax > prev_pos & imax < indxi); min_val]);

      coords(j,3) = 2*max_val - min_val - prev_min;
      prev_min = min_val;
      prev_pos = indxi;

      coords(j,4) = min_val;
      coords(j,5) = max_val;
    end

    all_coords = [all_coords; coords];
  end

  %orig_coords = all_coords;

  mval = nanmedian(all_coords(:,3:5), 1);
  sval = 1.4826 * mad(all_coords(:,3:5), 0, 1);

  if (any(isnan(mval) | isnan(sval)))
    all_coords = NaN(1,2);

    return;
  end

  all_coords = all_coords(all_coords(:,4) <= mval(2) + 3*sval(2), :);

  %rval = range(img(:));
  %keyboard

  if (ampl < 0)
    ampl = ((nanmax(all_coords(:,5)) - nanmin(all_coords(:,4))) - ...
           (mval(3)+sval(3) - (mval(2)-sval(2)))) / 15;
  end

  bads = (all_coords(:,3) < ampl);
  prev = [bads(2:end); false];

  all_coords(prev,:) = (all_coords(bads,:) + all_coords(prev,:))*0.5;
  all_coords(bads & ~prev, :) = [];
  %all_coords = all_coords(:,1:2);

  dists = bsxfun(@minus, all_coords(:,1), all_coords(:,1).').^2 + bsxfun(@minus, all_coords(:,2), all_coords(:,2).').^2;
  dists = triu(dists, 1);
  dists(dists == 0) = Inf;

  [bads, prev] = find(dists < (4*rads)^2);

  all_coords(prev,:) = (all_coords(bads,:) + all_coords(prev,:))*0.5;
  all_coords(bads(~ismember(bads, prev)), :) = [];
  all_coords = all_coords(:,1:2);

  sizes = size(all_coords, 1);

  %{
  figure;
  imagesc(orig_img);
  hold on;
  plot(systems(:,1), systems(:,2), 'w');
  plot(orig_sys(:,1), orig_sys(:,2), 'y');
  scatter(orig_coords(:,1), orig_coords(:,2), 'r');
  scatter(all_coords(:,1), all_coords(:,2), 'k');

  keyboard
  %}

  return;
end
