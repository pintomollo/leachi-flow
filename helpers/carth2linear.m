function [linear, total_distance] = carth2linear(pts_x, pts_y, circular)

  linear = [];
  total_distance = NaN;

  if (isempty(pts_x))
    return;
  end

  if (nargin < 2)
    circular = false;
    pts_y = pts_x(:,2);
    pts_x = pts_x(:,1);
  elseif (nargin < 3)
    if (islogical(pts_y))
      circular = pts_y;
      pts_y = pts_x(:,2);
      pts_x = pts_x(:,1);
    else
      circular = false;
    end
  end

  npts = numel(pts_x);

  if (circular && (pts_x(1) ~= pts_x(end) | pts_y(1) ~= pts_y(end)))
    pts_x = [pts_x; pts_x(1)];
    pts_y = [pts_y; pts_y(1)];
  end

  distance = sqrt(diff(pts_x).^2 + diff(pts_y).^2);
  cum_dist = [0; cumsum(distance)];

  total_distance = cum_dist(end);

  if (~(total_distance > 0))
    return;
  end

  if (nargout > 1)
    linear = cum_dist / total_distance;
  else
    linear = cum_dist;
  end

  if (numel(linear) > npts)
    linear = linear(1:npts);
  end

  return;
end
