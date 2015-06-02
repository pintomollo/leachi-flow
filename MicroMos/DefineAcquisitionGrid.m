function grid = DefineAcquisitionGrid(parameters)

  numberCorners = 150;
  flag_Harris = 0;
  thresh = 0.5;
  if (flag_Harris)
    method = 'Harris';
  else
    method = 'MinimumEigenvalue';
  end

  nelems = length(parameters.ImageIndexs);
  sizes = factorize_grid(nelems);

  nperms = size(sizes, 1);

  tmp_indx = NaN(nperms, 1);
  tmp_indx(1:2:end) = [1:ceil(nperms/2)];
  tmp_indx(2:2:end) = [nperms:-1:ceil(nperms/2)+1];

  sizes = sizes(tmp_indx, :);

  if (nelems == 1)
    trivial
  end

  props = NaN(nelems, 3);

  for i=1:nelems
    edges = imadm(load_img(i), false);

    [mval, stdval] = mymean(edges(:));
    props(i, 1) = mval + stdval;
  end

  img_size = size(edges);
  indexes = [1:nelems];

  best = 1;
  for i=size(sizes, 1):-1:2
    curr_indx = reshape(indexes, sizes(i,[2 1]));

    if (parameters.GridMode==0)
      curr_indx = permute(curr_indx, [2 1]);
    end

    data = props(curr_indx);

    vert_dist = data(1:end-1,:) + data(2:end, :);
    [rindx, cindx] = find(vert_dist==max(vert_dist(:)), 1);
    vindx = [curr_indx(rindx, cindx) curr_indx(rindx+1, cindx)];

    img1 = load_img(vindx(1));
    img2 = load_img(vindx(2));
    [v1, v2] = correl_edges(img1.', img2.', thresh);
    overlap = max(v1, v2);
    vpoints = corner(img1(end-overlap:end, :), method, numberCorners);
    vpoints(:,2) = vpoints(:,2) + (img_size(1)-overlap);

    horz_dist = data(:,1:end-1) + data(:,2:end);
    [rindx, cindx] = find(horz_dist==max(horz_dist(:)), 1);
    hindx = [curr_indx(rindx, cindx) curr_indx(rindx, cindx+1)];

    img1 = load_img(hindx(1));
    img2 = load_img(hindx(2));
    [h1, h2] = correl_edges(img1, img2, thresh);
    overlap = max(h1, h2);
    hpoints = corner(img1(:,end-overlap:end), method, numberCorners);
    hpoints(:,1) = hpoints(:,1) + (img_size(2)-overlap);

    keyboard
  end

  return;

  function [img] = load_img(img_indx)
    strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(img_indx));
    if strcmp(parameters.ImageFormat, '.mat')
        img = load(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum parameters.ImageFormat]));
        img = cell2mat(struct2cell(img));
    else
        img = imread(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum parameters.ImageFormat]));
    end
    img = img(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);

    if (size(img, 3) > 1)
      img = rgb2gray(img);
    end

    img = double(img);
  end
end

function [max1, max2] = correl_edges(img1, img2, thresh)

  img1 = bsxfun(@minus, img1, mean(img1, 1));
  img2 = bsxfun(@minus, img2, mean(img2, 1));

  corr1 = sum(bsxfun(@times, img1(:,end), img2), 1) ./ sqrt(sum(img1(:,end).^2, 1) * sum(img2.^2, 1));
  corr2 = sum(bsxfun(@times, img1, img2(:,1)), 1) ./ sqrt(sum(img1.^2, 1) * sum(img2(:,1).^2, 1));

  [val, max1] = max(corr1);

  if (val < thresh)
    max1 = NaN;
  end

  [val, max2] = max(corr2);

  if (val < thresh)
    max2 = NaN;
  else
    max2 = length(corr2) - max2 + 1;
  end

  return;
end
