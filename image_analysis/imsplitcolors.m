function img = imsplitcolors(img, dist)

  if (nargin == 1)
    dist = 5;
  end

  nbins = 255;

  if (~isa(img, 'uint8'))
    img = uint8(img*255);
  end

  [h,w,c] = size(img);

  img = rgb_hsv_mex(img);

  n2 = imsplitcolors_mex(img, zeros(nbins, 1));
  n2 = [n2; n2];

  n2 = colfilt(n2, [2*dist+1 1], 'sliding', @(y)(mean(y, 1)));

  [xmax, imax] = local_extrema(n2, dist);

  goods = (imax > dist & imax <= length(n2) - dist);
  xmax = xmax(goods);
  imax = imax(goods);

  [imax, indxs] = unique(mod(imax, nbins));
  xmax = xmax(indxs);

  [xmax, indxs] = sort(xmax);
  imax = imax(indxs(end:-1:1));

  imax = (imax(1:min(3, end))-1)/nbins;

  img = imsplitcolors_mex(img, imax);

  img = rgb_hsv_mex(img);

  return;
end
