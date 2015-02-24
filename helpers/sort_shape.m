function shapes = sort_shape(pts)

  if (isempty(pts))
    shapes = {};
    return;
  end

  dist = bsxfun(@minus, pts(:,1), pts(:,1).').^2 + ...
         bsxfun(@minus, pts(:,2), pts(:,2).').^2;

  neigh = sum(dist <= 2, 2) - 1;
  pts = [pts neigh];

  cross = (neigh > 2);
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

  return;
end
