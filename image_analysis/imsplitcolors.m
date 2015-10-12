function [img, nhist] = imsplitcolors(img, imax, nhist)

  dist = 10;
  if (nargin == 1)
    imax = [];
    nhist = [];
  elseif (nargin == 2)
    nhist = [];
  end

  if (~isnumeric(imax) || ~isfinite(imax))
    imax = [];
  end

  if (numel(imax) == 1)
    dist = imax;
    imax = [];
  end

  if (size(img, 3) < 3)
    if (~isempty(nhist) || ~isempty(imax))
      warning('Input arguments are not consistent.');
      img = [];
      nhist = [];
      return;
    end

    nhist = img;
    img = [];
  end

  if (isempty(nhist))
    nbins = 255;
    nhist = zeros(nbins, 1);
  else
    nbins = size(nhist);
  end

  compute_split = isempty(imax);

  if (~isempty(img))
    if (~isa(img, 'uint8'))
      img = uint8(img*255);
    end

    img = rgb_hsv_mex(img);

    if (compute_split)
      nhist = imsplitcolors_mex(img, [], nhist);
    end
  end

  if (nargout==1)
    if (compute_split)
      nhist = [nhist; nhist];

      ssize = size(nhist);
      ndim = length(ssize);

      if (ndim == 2 && (ssize(1) == 1 || ssize(2) == 1))
        ndim = 1;
        ssize = ssize(1);
      end

      dist = [dist(:).' ones(1, ndim-length(dist))*dist(1)];
      dist = dist(1:ndim);

      kernel = ones([2*dist+1 ones(1, 2-length(dist))]);
      kernel = kernel/numel(kernel);

      nhist = convn(nhist, kernel, 'same');
      %nhist = colfilt(nhist, [2*dist+1 1], 'sliding', @(y)(mean(y, 1)));

      [xmax, imax] = find_extrema(nhist, dist);

      goods = (imax(:,1) > dist(1)+1) & (imax(:,1) < ssize(1) - dist(1));
      %goods = (all(bsxfun(@gt, imax, dist+1), 2) & ...
      %         all(bsxfun(@lt, imax, ssize - dist), 2));

      if (ndim > 1)
        goods = goods & (imax(:,end) / ssize(2) < 0.91);
      end

      xmax = xmax(goods);
      imax = imax(goods, :);

      [imax, indxs] = unique(bsxfun(@mod, imax, nbins), 'rows');
      xmax = xmax(indxs);

      [xmax, indxs] = sort(xmax);
      imax = imax(indxs(end:-1:1),:);

      imax = bsxfun(@rdivide, imax(1:min(3, end),:)-1, nbins);
      imax = [imax; bsxfun(@times, ones(3-size(imax,1), ndim), imax(1,:))];
    end

    if (~isempty(img))
      img = imsplitcolors_mex(img, imax);

      img = rgb_hsv_mex(img);
    else
      img = imax;
    end
  end

  return;
end
