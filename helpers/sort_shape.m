function centers = sort_shape(pts, min_branch)

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
  sorted = NaN(2*npts, npos);

  gap = true;
  skip = -1;
  for i=1:npts
    if (gap)
      indx = find(tips, 1);
      val = 0;

      if (isempty(indx))
        [junk, indx] = min(dist(prev_indx, :), [], 2);
        [junk, indx_link] = min(cross(indx, :), [], 2);
        indx = cross_indx(indx_link);
      end
    else
      [val, indx] = min(dist(prev_indx, :), [], 2);
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
    else
      if (val > 2)
        sorted(i+skip, :) = NaN;
        skip = skip - 1;
        indx = prev_indx;
        link = pts(indx, :);
        gap = true;
      elseif (link(end)~=2)
        gap = true;
      end
    end
    sorted(i+skip, :) = link;

    prev_indx = indx;

    dist(:, prev_indx) = Inf;
    tips(prev_indx) = false;
  end

  if (link(end)==2)
    skip = skip + 1;
    [val_link, indx_link] = min(cross(indx, :), [], 2);
    sorted(i+skip, :) = pts(cross_indx(indx_link), :);
  end

  sorted = sorted(1:npts+skip,:);
  gaps = find(isnan(sorted(:,1)));
  ngaps = length(gaps);

  shapes = cell(ngaps+1, 1);
  indx = 1;
  for i=1:ngaps
    shapes{i} = sorted(indx:gaps(i)-1,:);
    indx=gaps(i)+1;
  end
  shapes{end} = sorted(indx:end,:);

  nbranches = ngaps+1;
  all_cross(:,end) = NaN;

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

    %scatter(pts(:,1), pts(:,2), 'MarkerEdgeColor', colors(i,:));
    lsq_line = [pts(:,1) ones(size(pts,1), 1)] \ pts(:,2);
    x_range = [min(pts(:,1)) max(pts(:,1))];
    %plot(x_range, x_range*lsq_line(1) + lsq_line(2), 'Color', colors(i,:)*0.5);

    if (pts(1,end)>2)
      all_cross(all(bsxfun(@eq,all_cross(:,1:2),pts(1,1:2)),2), end) = i;
      all_cross(end+1,:) = [pts(1,1) pts(1,1)*lsq_line(1) + lsq_line(2) i];
    end
    if (pts(end,end)>2)
      all_cross(all(bsxfun(@eq,all_cross(:,1:2),pts(end,1:2)),2), end) = i;
      all_cross(end+1,:) = [pts(end,1) pts(end,1)*lsq_line(1) + lsq_line(2) i];
    end

    all_lines(i,:) = [lsq_line.' x_range];
  end

  dist = bsxfun(@minus, all_cross(:,1), all_cross(:,1).').^2 + ...
         bsxfun(@minus, all_cross(:,2), all_cross(:,2).').^2;

  ncross = size(all_cross, 1);
  for i=1:ncross
    goods = (dist(:,i) <= 8);

    if (any(goods))

      links = all_cross(goods, [1 end]);
      goods = goods | ismember(all_cross(:, [1 end]), links, 'rows');

      new_center = mean(all_cross(goods,1:2), 1);
      all_cross(goods,1) = new_center(1);
      all_cross(goods,2) = new_center(2);

      dist(:, goods) = Inf;
    end
  end

  all_cross = unique(all_cross(~any(isnan(all_cross),2),:), 'rows');

  centers = NaN(0,2);
  for i=1:nbranches
    pts = shapes{i};

    if (~isempty(pts))
      tips = (all_cross(:,end)==i);
      switch sum(tips)
        case 0
          centers = [centers; all_lines(i,3)*[0 all_lines(i,1)] + all_lines(i,2); ...
                              all_lines(i,4)*[0 all_lines(i,1)] + all_lines(i,2)];
        case 1
          center = all_cross(tips,1:2);
          pts = bsxfun(@minus, pts(:,1:2), center);
          lsq_line = pts(:,1) \ pts(:,2);

          x_range = [min(pts(:,1)) max(pts(:,1))];
          [junk, indx] = max(abs(x_range));
          x_range = x_range(indx);

          centers = [centers; center; ([x_range lsq_line*x_range]) + center];
        case 2
          centers = [centers; all_cross(tips,1:2)];
        otherwise
          error('why?')
      end
      centers = [centers; NaN(1,2)];
    end
  end
  centers = centers(1:end-1,:);

  centers = polarize_flow(centers);

  return;
end
