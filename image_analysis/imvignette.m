function [estim] = imvignette(img, stack)
% IMVIGNETTE computes the vignette function according to Piccinini et al.

  if (nargin > 1)
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
    bkg(bw) = 0;

    if (isempty(stack))
      stack = cat(3, bkg, double(~bw));
    else
      stack(:,:,1) = stack(:,:,1) + bkg;
      stack(:,:,2) = stack(:,:,2) + double(~bw);
    end

    estim = stack;
  else

    estim = img(:,:,1) ./ img(:,:,2);
    estim = inpaint_nans(estim);

    %% Creates artifacts at the edges, where its most important
    %estim = gaussian_mex(estim, 15);
  end

  return;
end
