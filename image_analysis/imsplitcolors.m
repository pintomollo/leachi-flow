function imsplitcolors(img)

  dist = 5;
  nbins = 255;

  [h,w,c] = size(img);

  rgb = reshape(img, [h*w c]);
  hsv = rgb2hsv(rgb(:,1), rgb(:,2), rgb(:,3));

  cols = [0:nbins]/nbins;
  cols(end) = cols(end) + eps;

  [n, bins] = histc(hsv(:,1), cols);

  for i=1:nbins+1
    n(i) = sum(hsv(bins==i,2) .* hsv(bins==i,3));
  end

  n2 = [n(1:end-1); n(1:end-1)];
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

  dists = bsxfun(@minus, hsv(:,1), imax(:).');
  dists = cat(3, dists, dists+1, dists-1);
  [junk, indxs] = min(abs(dists), [], 3);

  mapping = dists(:,:,1);
  for i=2:3
    vals = dists(:,:,i);
    mapping(indxs==i) = vals(indxs==i);
  end

  objs = [0 1/3 2/3].';
  [junk, indxs] = min(abs(mapping), [], 2);

  mapping = mapping([1:size(mapping,1)].' +size(mapping,1)*(indxs-1));
  mapping = mapping + objs(indxs);
  mapping = mod(mapping, 1);

  hsv2 = hsv;
  hsv2(:,1) = mapping;
  rgb2 = hsv2rgb(hsv2);

  n3 = histc(hsv2(:,1), cols);
  n4 = [n3(1:end-1); n3(1:end-1)];

  figure;
  subplot(3,2,1);imagesc(img);
  subplot(3,2,2);imagesc(reshape(rgb2, [h w c]));
  subplot(3,2,3);plot(n2);hold on;scatter(imax*nbins + 1, n2(imax*nbins + 1), 'r');
  subplot(3,2,4);plot(n4);
  subplot(3,2,5);imagesc(rgb2gray(img))
  subplot(3,2,6);imagesc(rgb2gray(reshape(rgb2, [h w c])));

  return;
end
