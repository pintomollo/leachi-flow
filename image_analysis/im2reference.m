function [img, mask, hdil] = im2reference(img, minsize, hdil)

  if (nargin<2)
    img = size(img);
    minsize = 50;
    hdil = [];
  elseif (nargin < 3)
    img = size(img);
    hdil = [];
  end

  if (size(img,1) == 1)
    minrad = ceil(max(img)/minsize);
    minsize = (minrad^2);
    hdil = strel('disk', ceil(minrad/20));

    img = minrad;
    mask = minsize;
  else

    if (size(img,3) > 1)
      img = rgb2gray(img);
    end

    img = adjust_tif(img);
    th = graythresh(img(:))*max(img(:));
    mask = (img > th);
    mask = imclose(mask, hdil);
    mask = bwareaopen(mask, minsize);
    mask = imdilate(mask, hdil);

    if (any(mask(:)))
      img(~mask) = 0;
    end
    img = imnorm(img);
  end

  return;
end
