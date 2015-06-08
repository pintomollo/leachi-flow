function [xShift, yShift] = ShiftByCornerClustering(I1, I2, method, numberCorners, threshold)

  % Maybe we got a list of corners already
  if (size(I1, 2)==2)
    corners1 = I1;
  else
    tmp_img = double(I1);
    tmp_img(isnan(I1)) = nanmean(tmp_img(:));
    corners1 = corner(tmp_img, method, numberCorners);
  end

  % Maybe we got a list of corners already
  if (size(I2, 2)==2)
    corners2 = I2;
  else
    tmp_img = double(I2);
    tmp_img(isnan(I2)) = nanmean(tmp_img(:));
    corners2 = corner(tmp_img, method, numberCorners);
  end

  % Compute the all-to-all distances
  dx = bsxfun(@minus, corners1(:,1), corners2(:,1).');
  dy = bsxfun(@minus, corners1(:,2), corners2(:,2).');

  % Cluster them
  [clusts, nums] = cluster_vector_mex(dx(:), dy(:), threshold);

  % Get the biggest cluster
  [vals, indxs] = sort(nums, 'descend');

  % Average them to get the corresponding shift
  xShift = mean(dx(clusts==indxs(1)));
  yShift = mean(dy(clusts==indxs(1)));

  return;
end
