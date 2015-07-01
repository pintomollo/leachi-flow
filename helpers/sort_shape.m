function centers = sort_shape(pts, min_branch, min_curve)

  if (isempty(pts))
    centers = NaN(0,2);
    return;
  end

  if (nargin < 2)
    min_branch = 2;
  end

  dist = bsxfun(@minus, pts(:,1), pts(:,1).').^2 + ...
         bsxfun(@minus, pts(:,2), pts(:,2).').^2;

  neigh = sum(dist <= 2, 2) - 1;
  pts = [pts neigh];

  cross = (neigh > 2);
  all_cross = pts(cross,:);
  cross_indx = find(cross);
  cross = dist(:,cross);

  tips = (neigh == 1);
  tips_indx = find(tips);

  [npts, npos] = size(pts);
  npts = npts + sum(all_cross(:,end)-1);
  sorted = NaN(2*npts, npos);

  %figure;hold on;
  %scatter(pts(:,1), pts(:,2), 'b')
  %scatter(all_cross(:,1), all_cross(:,2), 'r');
  %scatter(pts(tips,1), pts(tips,2), 'k');

  gap = true;
  skip = -1;
  for i=1:npts
    if (gap)
      indx = find(tips, 1);
      val = 0;

      if (isempty(indx))
        [junk, indx] = min(dist(prev_indx, :), [], 2);

        if (~isfinite(junk))
          break;
        end

        [junk, indx_link] = min(cross(indx, :), [], 2);
        indx = cross_indx(indx_link);
      end
    else
      [val, indx] = min(dist(prev_indx, :), [], 2);
      if (val>2)
        [junk, indx_link] = min(cross(prev_indx, :), [], 2);
        indx = cross_indx(indx_link);
      end
    end
    link = pts(indx, :);

    if (gap)
      gap = false;
      skip = skip + 1;

      if (link(end)>1)
        [val_link, indx_link] = min(cross(indx, :), [], 2);
        sorted(i+skip, :) = pts(cross_indx(indx_link), :);
        skip = skip + 1;
      end
    elseif (link(end)~=2)
      gap = true;
    end
    sorted(i+skip, :) = link;

    prev_indx = indx;

    dist(:, prev_indx) = Inf;
    tips(prev_indx) = false;

    %scatter(link(1), link(2), 'm');
    %drawnow
  end

  if (link(end)==2 && ~gap)
    skip = skip + 1;
    [val_link, indx_link] = min(cross(indx, :), [], 2);
    sorted(i+skip, :) = pts(cross_indx(indx_link), :);
  end

  last = find(~isnan(sorted(:,1)), 1, 'last');
  sorted = sorted(1:last,:);
  gaps = find(isnan(sorted(:,1)));
  ngaps = length(gaps);

  dist = bsxfun(@minus, all_cross(:,1), all_cross(:,1).').^2 + bsxfun(@minus, all_cross(:,2), all_cross(:,2).').^2;
  ncross = size(dist, 1);
  fuse = eye(ncross);
  for i=1:ncross
    goods = (dist(i,:) <= 2);
    grp = goods | fuse(i,:);
    fuse(grp, grp) = true;
  end
  fuse = logical(unique(fuse, 'rows'));

  for i=1:size(fuse,1)
    vals = all_cross(fuse(i,:),:);
    replace = ismember(sorted, vals, 'rows');
    news = mean(vals,1);
    sorted(replace,1) = news(1);
    sorted(replace,2) = news(2);
  end
  all_cross = unique(sorted(sorted(:,3)>2,1:2), 'rows');

  %colors = jet(ngaps+1);
  %figure;
  %hold on;

  shapes = cell(ngaps+1, 1);
  indx = 1;
  for i=1:ngaps
    shapes{i} = sorted(indx:gaps(i)-1,:);
    indx=gaps(i)+1;
    %scatter(shapes{i}(:,1),shapes{i}(:,2), 'MarkerEdgeColor', colors(i,:));
  end
  shapes{end} = sorted(indx:end,:);
  %scatter(shapes{end}(:,1),shapes{end}(:,2), 'MarkerEdgeColor', colors(end,:));
  %plot(sorted(:,1), sorted(:,2), 'k');

  tmp_shapes = {};
  for i=1:length(shapes)
    pts = shapes{i};
    pacs = impac(pts, min_curve, 'recursive');

    knots = find(ismember(pts, pacs, 'rows'));
    pts(knots(2:end-1), 3) = 3;
    for j=1:length(knots)-1
      tmp_shapes{end+1} = pts(knots(j):knots(j+1),:);
    end
  end
  shapes = tmp_shapes;
  all_pts = cat(1, shapes{:});

  nbranches = length(shapes);
  orig_cross = unique(all_pts(all_pts(:,3)>2,1:2), 'rows');
  all_cross = NaN(size(orig_cross) + [0 nbranches]);
  all_cross(:,1:2) = orig_cross;

  %scatter(orig_cross(:,1), orig_cross(:,2), 'k');

  %colors = jet(nbranches);
  %figure;
  %hold on;

  all_lines = NaN(nbranches, 4);
  for i=1:nbranches
    pts = shapes{i};

    if (size(unique(pts,'rows'),1) < min_branch)
      shapes{i} = [];
      continue;
    end

    ranges = myregress(pts);
    %scatter(pts(:,1), pts(:,2), 'MarkerEdgeColor', colors(i,:));
    %plot(ranges([1 3]), ranges([2 4]), 'Color', colors(i,:)*0.5);

    if (pts(1,end)>2)
      cindx = all(bsxfun(@eq,orig_cross,pts(1,1:2)),2);
      rindx = find(isnan(all_cross(cindx,:)), 1, 'first');
      all_cross(cindx, rindx) = i;
      all_cross(cindx, 1:2) = all_cross(cindx, 1:2) + ranges(1:2);
    end
    if (pts(end,end)>2)
      cindx = all(bsxfun(@eq,orig_cross,pts(end,1:2)),2);
      rindx = find(isnan(all_cross(cindx,:)), 1, 'first');
      all_cross(cindx, rindx) = i;
      all_cross(cindx, 1:2) = all_cross(cindx, 1:2) + ranges(3:4);
    end

    all_lines(i,:) = ranges;
  end
  ncross = sum(~isnan(all_cross(:,3:end)), 2)+1;
  all_cross(:,1:2) = bsxfun(@rdivide, all_cross(:,1:2), ncross);
  %scatter(all_cross(:,1), all_cross(:,2), 'k');

  %colors = jet(nbranches);
  %figure;
  %hold on;

  centers = NaN(0,2);
  for i=1:nbranches
    pts = shapes{i};

    if (~isempty(pts))
      tips = (any(all_cross(:,3:end)==i, 2));
      switch sum(tips)
        case 0
          centers = [centers; all_lines(i,1:2); all_lines(i,3:4)];
        case 1
          center = all_cross(tips,1:2);
          ranges = myregress(pts, center);
          centers = [centers; ranges(1:2); ranges(3:4)];
        case 2
          centers = [centers; all_cross(tips,1:2)];
        otherwise
          error('why?')
      end
      centers = [centers; NaN(1,2)];

      %scatter(pts(:,1), pts(:,2), 'MarkerEdgeColor', colors(i,:));
      %plot(centers(end-2:end-1,1), centers(end-2:end-1,2), 'Color', colors(i,:)*0.5);
    end
  end
  centers = centers(1:end-1,:);
  %scatter(all_cross(:,1), all_cross(:,2), 'k');

  %keyboard
  centers = polarize_flow(centers);

  return;
end

function extrema = myregress(pts, pivot)

  if (nargin < 2)
    ranges = [min(pts, [], 1); max(pts, [], 1)];
    rval = diff(ranges, 1, 1);

    if (rval(1) > rval(2))
      lsq_line = [pts(:,1) ones(size(pts,1), 1)] \ pts(:,2);
      extrema = [ranges(:,1).'; (ranges(:,1).')*lsq_line(1) + lsq_line(2)];
      if (pts(1,1)>pts(end,1))
        extrema = extrema(:,[2 1]);
      end
    else
      lsq_line = [pts(:,2) ones(size(pts,1), 1)] \ pts(:,1);
      extrema = [(ranges(:,2).')*lsq_line(1) + lsq_line(2); ranges(:,2).'];
      if (pts(1,2)>pts(end,2))
        extrema = extrema(:,[2 1]);
      end
    end
  else
    pts = bsxfun(@minus, pts(:,1:2), pivot);
    ranges = [min(pts, [], 1); max(pts, [], 1)];
    [junk, indxs] = min(abs(ranges), [], 1);
    ranges(indxs(1), 1) = 0;
    ranges(indxs(2), 2) = 0;

    rval = diff(ranges, 1, 1);

    if (rval(1) > rval(2))
      lsq_line = pts(:,1) \ pts(:,2);
      extrema = [ranges(:,1).' + pivot(1); (ranges(:,1).')*lsq_line(1) + pivot(2)];
      if (pts(1,1)>pts(end,1))
        extrema = extrema(:,[2 1]);
      end
    else
      lsq_line = pts(:,2) \ pts(:,1);
      extrema = [(ranges(:,2).')*lsq_line(1) + pivot(1); ranges(:,2).' + pivot(2)];
      if (pts(1,2)>pts(end,2))
        extrema = extrema(:,[2 1]);
      end
    end
  end
  extrema = extrema(:).';

  return;
end
