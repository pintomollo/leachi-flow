function [estim, stack] = imvignette(img, stack)

  if (size(img, 3) > 1)
    img = rgb2gray(img);
  end

  img = double(img);
  edges = imadm(img, false);
  thresh = graythresh(edges);

  bw = (edges > 0.5*thresh);
  bw = imdilate(bw, strel('disk', 21));
  bw = imerode(bw, strel('disk', 9));
  bw = bwareaopen(bw, 225);

  bkg = img;
  bkg(bw) = NaN;

  if (nargin < 2 || isempty(stack))
    stack = bkg;
  else
    stack = cat(3, stack, bkg);
  end

  estim = nanmean(stack, 3);
  estim = inpaint_nans(estim);
  estim = gaussian_mex(estim, 15);

  return;
end
