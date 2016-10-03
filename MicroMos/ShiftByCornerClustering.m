function [xShift, yShift, nums] = ShiftByCornerClustering(I1, I2, method, numberCorners, threshold)

  dx = [];

  % Maybe we got the distances straight
  if (size(I1, 2) == 1)
    dx = I1;
  % Maybe we got a list of corners already
  elseif (size(I1, 2)==2)
    corners1 = I1;
  else
    tmp_img = double(I1);
    tmp_img(isnan(I1)) = nanmean(tmp_img(:));
    corners1 = GriddedCorner(tmp_img, method, numberCorners);
  end

  % Maybe we got the distances straight
  if (size(I2, 2) == 1)
    dy = I2;
  % Maybe we got a list of corners already
  elseif (size(I2, 2)==2)
    corners2 = I2;
  else
    tmp_img = double(I2);
    tmp_img(isnan(I2)) = nanmean(tmp_img(:));
    corners2 = GriddedCorner(tmp_img, method, numberCorners);
  end

  % Compute the all-to-all distances if need be
  if (isempty(dx))
    dx = bsxfun(@minus, corners1(:,1), corners2(:,1).');
    dy = bsxfun(@minus, corners1(:,2), corners2(:,2).');
  end

  % Cluster them
  [clusts, nums] = cluster_vector_mex(dx(:), dy(:), threshold);

  % Get the biggest cluster
  [vals, indxs] = sort(nums, 'descend');

  % Average them to get the corresponding shift
  if (nargout == 2)
    xShift = mean(dx(clusts==indxs(1)));
    yShift = mean(dy(clusts==indxs(1)));
  else
    xShift = mymean(dx,1,clusts);
    yShift = mymean(dy,1,clusts);

    xShift = xShift(indxs);
    yShift = yShift(indxs);
    nums = vals;
  end

  return;
end
