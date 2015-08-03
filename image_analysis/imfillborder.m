function [img, bkgs] = imfillborder(img)

  [h,w,c] = size(img);

  bkgs = NaN(1, c);

  borders = (any(isnan(img), 3) | all(img==0, 3));
  borders = ~(imfill(~borders, 'holes'));

  empty_rows = all(borders, 2);
  empty_cols = all(borders, 1);
  real_h = sum(~empty_rows);
  real_w = sum(~empty_cols);

  first_h = find(~empty_rows, 1, 'first');
  last_h = find(~empty_rows, 1, 'last');
  first_w = find(~empty_cols, 1, 'first');
  last_w = find(~empty_cols, 1, 'last');

  rim = ceil([real_h, real_w]/100);
  rim_image = true([h, w]);
  rim_image(rim(2)+first_h:last_h-rim(2)+1,rim(1)+first_w:last_w-rim(1)+1) = false;

  for i = 1:c
    tmp_img = img(:,:,i);
    bkgs(i) = nanmedian(double(tmp_img(rim_image & ~borders)));
    tmp_img(borders) = bkgs(i);
    img(:,:,i) = tmp_img;
  end

  if (nargout > 1)
    bkgs = cast(bkgs, class(img));
  end

  return;
end
