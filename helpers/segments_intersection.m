function [valids, pos] = segments_intersection(segm, tests)

  xs = [tests(1:end-1,1) tests(2:end, 1)];
  ys = [tests(1:end-1,2) tests(2:end, 2)];

  origin = segm(1,:);

  segm = segm(2,:) - origin;
  len = sqrt(sum(segm.^2));

  tmp_x = xs - origin(1);
  tmp_y = ys - origin(2);

  angles = segm / len;

  rot_x = tmp_x*angles(1) + tmp_y*angles(2);
  rot_y = tmp_y*angles(1) - tmp_x*angles(2);

  dists = rot_x(:,2) + (rot_x(:,1) - rot_x(:,2)) .* rot_y(:,2) ./ (diff(rot_y, [], 2));
  valids = (dists > 0 & dists < len & xor(rot_y(:,1) > 0, rot_y(:,2) > 0));

  pos = NaN(length(valids), 2);

  if (any(valids))
    pos(valids, :) = bsxfun(@plus, bsxfun(@times, dists(valids), angles), origin);
  end

  return;
end
