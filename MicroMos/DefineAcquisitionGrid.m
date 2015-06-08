function [MatricesGLOBAL, images_ordering] = DefineAcquisitionGrid(parameters)

  numberCorners = 150;
  flag_Harris = 0;
  thresh = 0.5;
  pix_thresh = 5;
  if (flag_Harris)
    method = 'Harris';
  else
    method = 'MinimumEigenvalue';
  end

  MatricesGLOBAL = eye(3);
  images_ordering = [];

  nelems = length(parameters.ImageIndexs);
  if (nelems == 1)
    images_ordering = 1;

    return;
  end

  sizes = factorize_grid(nelems);

  nperms = size(sizes, 1);

  tmp_indx = NaN(nperms, 1);
  tmp_indx(1:2:end) = [1:ceil(nperms/2)];
  tmp_indx(2:2:end) = [nperms:-1:ceil(nperms/2)+1];

  sizes = sizes(tmp_indx, :);

  props = NaN(nelems, 1);
  all_edges = cell(nelems, 2);

  for i=1:nelems
    edges = imadm(load_img(i), false);
    all_edges{i,1} = nanmean(edges, 2).';
    all_edges{i,2} = nanmean(edges, 1);

    [mval, stdval] = mymean(edges(:));
    props(i, 1) = mval + stdval;
  end

  avg_edges = mean(props(:));
  img_size = size(edges);
  indexes = [1:nelems];

  all_shifts = NaN(size(sizes, 1), 6);

  best = 0;
  for i=size(sizes, 1):-1:1
    curr_indx = reshape(indexes, sizes(i,[2 1]));

    if (parameters.GridMode==0)
      curr_indx = permute(curr_indx, [2 1]);
    end

    data = props(curr_indx);
    if (any(sizes(i,:)==1))
      data = reshape(data, sizes(i,:));
    end

    if (sizes(i,1) > 1)
      vert_dist = data(1:end-1,:) + data(2:end, :);
      [rindx, cindx] = find(vert_dist==max(vert_dist(:)), 1);
      vindx = [curr_indx(rindx, cindx) curr_indx(rindx+1, cindx)];

      img1 = load_img(vindx(1));
      img2 = load_img(vindx(2));

      strong1 = all_edges{vindx(1), 1} > avg_edges;
      strong2 = all_edges{vindx(2), 1} > avg_edges;

      [estims, corr1] = correl_edges(img1.', img2.', strong1, strong2, thresh);
      overlap = round(mean(estims));
      if (isnan(overlap) || overlap < pix_thresh || overlap > img_size(1)-pix_thresh)
        continue
      end

      vpoints = corner(img1(end-overlap:end, :), method, numberCorners);
      vpoints(:,2) = vpoints(:,2) + (img_size(1)-overlap);

      vpoints2 = corner(img2(1:overlap, :), method, numberCorners);
      [vshiftx, vshifty] = ShiftByCornerClustering(vpoints, vpoints2, method, numberCorners, 1);

      vpoints2 = LKTracker(img1, img2, vpoints, [vshiftx vshifty]);
      indices2 = CheckPointsAndNAN(img2, vpoints2);
      vpoints2 = vpoints2(indices2,:);
      vpoints = vpoints(indices2,:);

      vshift = mean(vpoints - vpoints2);
      all_shifts(i,1:2) = vshift;
      all_shifts(i,5) = corr1;
    end

    if (sizes(i,2) > 1)
      horz_dist = data(:,1:end-1) + data(:,2:end);
      [rindx, cindx] = find(horz_dist==max(horz_dist(:)), 1);
      hindx = [curr_indx(rindx, cindx) curr_indx(rindx, cindx+1)];

      img1 = load_img(hindx(1));
      img2 = load_img(hindx(2));

      strong1 = all_edges{hindx(1), 2} > avg_edges;
      strong2 = all_edges{hindx(2), 2} > avg_edges;

      [estims, corr2] = correl_edges(img1, img2, strong1, strong2, thresh);
      overlap = round(mean(estims));

      if (isnan(overlap) || overlap < pix_thresh || overlap > img_size(2)-pix_thresh)
        continue
      end

      hpoints = corner(img1(:,end-overlap:end), method, numberCorners);
      hpoints(:,1) = hpoints(:,1) + (img_size(2)-overlap);

      hpoints2 = corner(img2(:,1:overlap), method, numberCorners);
      [hshiftx, hshifty] = ShiftByCornerClustering(hpoints, hpoints2, method, numberCorners, 1);

      hpoints2 = LKTracker(img1, img2, hpoints, [hshiftx hshifty]);
      indices2 = CheckPointsAndNAN(img2, hpoints2);
      hpoints2 = hpoints2(indices2,:);
      hpoints = hpoints(indices2,:);

      hshift = mean(hpoints - hpoints2);
      all_shifts(i,3:4) = hshift;
      all_shifts(i,6) = corr2;
    end

    best = i;
    break;
  end

  if (best > 0)
    best_size = sizes(best,:);

    gridhx = [0:best_size(2)-1]*all_shifts(best,3);
    gridhy = [0:best_size(2)-1]*all_shifts(best,4);

    gridvx = [0:best_size(1)-1]*all_shifts(best,1);
    gridvy = [0:best_size(1)-1]*all_shifts(best,2);

    gridx = bsxfun(@plus, gridhx, gridvx.');
    gridy = bsxfun(@plus, gridhy, gridvy.');

    curr_indx = reshape(indexes, best_size([2 1]));

    if (parameters.GridMode==0)
      curr_indx = permute(curr_indx, [2 1]);
    end

    MatricesGLOBAL = repmat(eye(3), [1 1 nelems]);
    MatricesGLOBAL(1,3, curr_indx(:)) = gridx(:);
    MatricesGLOBAL(2,3, curr_indx(:)) = gridy(:);

    curr_flip = flipud(curr_indx);
    curr_indx(:,2:2:end) = curr_flip(:,2:2:end);

    images_ordering = curr_indx(:);
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

function [maxs, max_corr] = correl_edges(img1, img2, goods1, goods2, thresh)

  maxs = [];

  img1 = bsxfun(@minus, img1, mean(img1, 1));
  img2 = bsxfun(@minus, img2, mean(img2, 1));

  corr1 = sum(bsxfun(@times, img1(:,end), img2), 1) ./ sqrt(sum(img1(:,end).^2, 1) * sum(img2.^2, 1));
  corr2 = sum(bsxfun(@times, img1, img2(:,1)), 1) ./ sqrt(sum(img1.^2, 1) * sum(img2(:,1).^2, 1));

  [val1, max1] = max(corr1);

  if (val1 >= thresh && mean(corr1 > thresh) <= 0.25)
    maxs = max1;
  elseif (any(goods1(1:end-1)))
    good_column = find(goods1(1:end-1), 1, 'last');

    corr1 = sum(bsxfun(@times, img1(:,good_column), img2), 1) ./ sqrt(sum(img1(:,good_column).^2, 1) * sum(img2.^2, 1));
    [val1, max1] = max(corr1);

    if (val1 >= thresh && mean(corr1 > thresh) <= 0.25)
      maxs = (length(corr1)-good_column) + max1;
    end
  end

  [val2, max2] = max(corr2);

  if (val2 >= thresh && mean(corr2 > thresh) <= 0.25)
    maxs = [maxs length(corr2) - max2 + 1];
  elseif (any(goods2(2:end)))
    good_column = find(goods2(2:end), 1, 'first')+1;

    corr2 = sum(bsxfun(@times, img1, img2(:,good_column)), 1) ./ sqrt(sum(img1.^2, 1) * sum(img2(:,good_column).^2, 1));

    [val2, max2] = max(corr2);

    if (val2 >= thresh && mean(corr2 > thresh) <= 0.25)
      maxs = [maxs length(corr2) - max2 + good_column];
    elseif (isempty(maxs))
      maxs = NaN;
    end
  end

  max_corr = max(val1, val2);

  return;
end
