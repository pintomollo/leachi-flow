function [img, nhist] = imsplitcolors(img, dist)

  nhist = [];
  imax = [];
  if (nargin == 1)
    dist = 5;
  elseif (numel(dist) > 1)
    if (numel(dist) > 3)
      nhist = dist;
      dist = 5;
    else
      imax = dist;
      dist = 5;
    end
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
    nhist = nhist(:);
    nbins = length(nhist);
  end

  compute_split = isempty(imax);

  if (~isempty(img))
    if (~isa(img, 'uint8'))
      img = uint8(img*255);
    end

    img = rgb_hsv_mex(img);

    if (compute_split)
      nhist = imsplitcolors_mex(img, nhist);
    end
  end

  if (nargout==1)
    if (compute_split)
      nhist = [nhist; nhist];

      nhist = colfilt(nhist, [2*dist+1 1], 'sliding', @(y)(mean(y, 1)));

      [xmax, imax] = local_extrema(nhist, dist);

      goods = (imax > dist & imax <= length(nhist) - dist);
      xmax = xmax(goods);
      imax = imax(goods);

      [imax, indxs] = unique(mod(imax, nbins));
      xmax = xmax(indxs);

      [xmax, indxs] = sort(xmax);
      imax = imax(indxs(end:-1:1));

      imax = (imax(1:min(3, end))-1)/nbins;
      imax = [imax ones(1, 3-length(imax))*imax(1)];
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
