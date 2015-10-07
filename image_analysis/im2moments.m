function [means, angle, stds] = im2moments(img, xx, yy, obj_angle)

  if (nargin < 2)
    [h,w,c] = size(img);

    xx = repmat(single([1:w]), h, 1);
    yy = repmat(single([1:h]).', 1, w);

    obj_angle = [];
  elseif (nargin < 3)
    [h,w,c] = size(img);

    obj_angle = xx(1);

    xx = repmat(single([1:w]), h, 1);
    yy = repmat(single([1:h]).', 1, w);
  elseif (nargin < 4)

    obj_angle = [];
  end

  if (size(img, 3) > 1)
    img = im2reference(img);
  end

  if (islogical(img))

    x = xx(img);
    y = yy(img);

    xbar = mean(x);
    ybar = mean(y);

    x = x - xbar;
    y = y - ybar;

    npts = length(x);

    % Calculate normalized second central moments for the region. 1/12 is
    % the normalized second central moment of a pixel with unit length.
    uxx = sum(x.^2)/npts + 1/12;
    uyy = sum(y.^2)/npts + 1/12;
    uxy = sum(x.*y)/npts;

    % Calculate orientation.
    if (uyy > uxx)
      num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
      den = 2*uxy;
    else
      num = 2*uxy;
      den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
    end
    if (num == 0) && (den == 0)
      angle = 0;
    else
      angle = -atan2(num, den);
    end

    if (~isempty(obj_angle))
      if (abs(angle - obj_angle) > pi/2)
        angle = angle + pi;
      end
    end

    if (nargout > 2)
      x2 = x*cos(angle) - y*sin(angle);
      y2 = x*sin(angle) + y*cos(angle);

      uxx = sum(x2.^2)/npts + 1/12;
      uyy = sum(y2.^2)/npts + 1/12;

      stds = double(sqrt([uxx uyy]));
    end

    means = double([xbar ybar]);
    angle = double(angle);
  else

    goods = (img>0);

    x = xx(goods);
    y = yy(goods);
    w = double(img(goods));
    w = single(w / sum(w));

    xbar = sum(x .* w);
    ybar = sum(y .* w);

    x = x - xbar;
    y = y - ybar;

    npts = length(x);

    % Calculate normalized second central moments for the region. 1/12 is
    % the normalized second central moment of a pixel with unit length.
    uxx = sum((x.^2) .* w) + 1/12;
    uyy = sum((y.^2) .* w) + 1/12;
    uxy = sum((x.*y) .* w);

    % Calculate orientation.
    if (uyy > uxx)
      num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
      den = 2*uxy;
    else
      num = 2*uxy;
      den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
    end
    if (num == 0) && (den == 0)
      angle = 0;
    else
      angle = -atan2(num, den);
    end

    if (~isempty(obj_angle))
      if (abs(angle - obj_angle) > pi/2)
        angle = angle + pi;
      end
    end

    if (nargout > 2)
      x2 = x*cos(angle) - y*sin(angle);
      y2 = x*sin(angle) + y*cos(angle);

      uxx = sum((x2.^2) .* w) + 1/12;
      uyy = sum((y2.^2) .* w) + 1/12;

      stds = double(sqrt([uxx uyy]));
    end

    means = double([xbar ybar]);
    angle = double(angle);
  end

  return;
end
