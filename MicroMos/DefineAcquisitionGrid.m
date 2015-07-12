function [MatricesGLOBAL, images_ordering] = DefineAcquisitionGrid(parameters)

  numberCorners = parameters.numberCorners;
  flag_Harris = 0;
  thresh = 0.3;
  pix_thresh = 5;
  prop_thresh = 0.25;
  if (flag_Harris)
    method = 'Harris';
  else
    method = 'MinimumEigenvalue';
  end

  MatricesGLOBAL = eye(3);
  images_ordering = [];

  nelems = length(parameters.ImageIndexs);
  if (nelems <= 1)
    images_ordering = 1;

    return;
  end

  sizes = factorize_grid(nelems);

  nperms = size(sizes, 1);

  tmp_indx = NaN(nperms, 1);
  tmp_indx(1:2:end) = [1:ceil(nperms/2)];
  tmp_indx(2:2:end) = [nperms:-1:ceil(nperms/2)+1];

  sizes = sizes(tmp_indx, :);
  sizes = sizes(~any(sizes==1,2),:);

  props = NaN(nelems, 4);
  all_edges = cell(nelems, 2);

  img_size = [];
  total_edges = 0;
  for i=1:nelems
    edges = imadm(load_img(i), false);
    if (isempty(img_size))
      img_size = size(edges);
      rim_size = ceil(img_size/4);
    end
    all_edges{i,1} = nanmean(edges, 2).';
    all_edges{i,2} = nanmean(edges, 1);

    props(i, 1) = mean(all_edges{i,1}(end-rim_size(1)+1:end));
    props(i, 2) = mean(all_edges{i,1}(1:rim_size(1)));
    props(i, 3) = mean(all_edges{i,2}(end-rim_size(2)+1:end));
    props(i, 4) = mean(all_edges{i,2}(1:rim_size(2)));

    total_edges = total_edges + mean(all_edges{i,1});
  end

  avg_edges = total_edges / nelems;
  indexes = [1:nelems];

  all_shifts = NaN(size(sizes, 1), 8);

  for i=size(sizes, 1):-1:1
    curr_indx = reshape(indexes, sizes(i,[2 1]));

    if (parameters.GridMode==0)
      curr_indx = permute(curr_indx, [2 1]);
    end

    data1 = props(curr_indx,1);
    data1 = reshape(data1, sizes(i,:));
    data2 = props(curr_indx,2);
    data2 = reshape(data2, sizes(i,:));

    if (sizes(i,1) > 1)
      vert_dist = data1(1:end-1,:) .* data2(2:end, :);

      [vert_edges, sub_ind] = sort(vert_dist(:), 'descend');
      [rindx, cindx] = ind2sub(size(vert_dist), sub_ind);

      ncands = length(rindx);

      corrs = zeros(1,img_size(1));
      for j=1:ncands
        if (vert_edges(j) < 3*(avg_edges^2) && j > 3)
          break;
        end

        vindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j)+1, cindx(j))];

        img1 = load_img(vindx(1));
        img2 = load_img(vindx(2));

        [img_corr] = correl_images(img1.', img2.', all_edges{vindx(1),1}, all_edges{vindx(2),1}, thresh, avg_edges);
        corrs = corrs + img_corr;
      end

      [corr1, overlap] = max(corrs);

      if (~any(corrs) || isnan(corr1) || overlap < pix_thresh || overlap > img_size(1)-pix_thresh)
        continue;
      end

      all_vpoints = cell([1 ncands]);
      all_vpoints2 = cell([1 ncands]);
      for j=1:ncands
        if (vert_edges(j) < 3*(avg_edges^2) && j > 3)
          break;
        end

        vindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j)+1, cindx(j))];

        img1 = load_img(vindx(1));
        img2 = load_img(vindx(2));

        vpoints = corner(img1(end-overlap+1:end, :), method, numberCorners);
        vpoints(:,2) = vpoints(:,2) + (img_size(1)-overlap);

        metric1 = cornermetric(img1, method);
        metric1 = metric1(sub2ind(img_size, vpoints(:,2), vpoints(:,1)));

        vpoints2 = corner(img2(1:overlap, :), method, numberCorners);
        metric2 = cornermetric(img2, method);
        metric2 = metric2(sub2ind(img_size, vpoints2(:,2), vpoints2(:,1)));

        all_vpoints{j} = [vpoints metric1 j*ones(size(metric1))];
        all_vpoints2{j} = [vpoints2 metric2 j*ones(size(metric2))];
      end

      tmp1 = cat(1, all_vpoints{:});
      tmp2 = cat(1, all_vpoints2{:});

      if (isempty(tmp1) || all(isnan(tmp1(:))))
        continue;
      end

      [tmp_v, tmp_indxs] = sort(tmp1(:,3), 'descend');
      tmp1 = tmp1(tmp_indxs(1:min(numberCorners, end)),:);

      [vshiftx, vshifty] = ShiftByCornerClustering(tmp1(:,1:2), tmp2(:,1:2), method, numberCorners, 1);

      for j=1:length(rindx)
        curr_pts = (tmp1(:,4) == j);

        if (any(curr_pts))

          vindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j)+1, cindx(j))];

          img1 = load_img(vindx(1));
          img2 = load_img(vindx(2));

          vpoints = tmp1(curr_pts, 1:2);

          vpoints2 = LKTracker(img1, img2, vpoints, [vshiftx vshifty]);
          indices2 = CheckPointsAndNAN(img2, vpoints2);
          vpoints2 = vpoints2(indices2,:);
          vpoints = vpoints(indices2,:);

          all_vpoints{j} = vpoints;
          all_vpoints2{j} = vpoints2;
        else
          all_vpoints{j} = NaN(0,2);
          all_vpoints2{j} = NaN(0,2);
        end
      end

      vpoints = cat(1, all_vpoints{:});
      vpoints2 = cat(1, all_vpoints2{:});

      if (isempty(vpoints) || all(isnan(vpoints(:))))
        continue;
      end

      vshift = mean(vpoints - vpoints2, 1);

      if (abs(vshift(1)/vshift(2)) > prop_thresh)
        continue;
      end

      all_shifts(i,1:2) = vshift;
      all_shifts(i,5) = corr1;
      all_shifts(i,7) = size(vpoints,1);
    else
      all_shifts(i,5) = 0;
      all_shifts(i,7) = 0;
    end

    data1 = props(curr_indx,3);
    data1 = reshape(data1, sizes(i,:));
    data2 = props(curr_indx,4);
    data2 = reshape(data2, sizes(i,:));

    if (sizes(i,2) > 1)
      horz_dist = data1(:,1:end-1) .* data2(:,2:end);

      [horz_edges, sub_ind] = sort(horz_dist(:), 'descend');
      [rindx, cindx] = ind2sub(size(horz_dist), sub_ind);

      ncands = length(rindx);

      corrs = zeros(1,img_size(2));
      for j=1:ncands
        if (horz_edges(j) < 3*(avg_edges^2) && j > 3)
          break;
        end

        hindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j), cindx(j)+1)];

        img1 = load_img(hindx(1));
        img2 = load_img(hindx(2));

        [img_corr] = correl_images(img1, img2, all_edges{hindx(1),2}, all_edges{hindx(2),2}, thresh, avg_edges);
        corrs = corrs + img_corr;
      end

      [corr2, overlap] = max(corrs);

      if (~any(corrs) || isnan(corr2) || overlap < pix_thresh || overlap > img_size(2)-pix_thresh)
        continue;
      end

      all_hpoints = cell([1 ncands]);
      all_hpoints2 = cell([1 ncands]);
      for j=1:ncands
        if (horz_edges(j) < 3*(avg_edges^2) && j > 3)
          break;
        end

        hindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j), cindx(j)+1)];

        img1 = load_img(hindx(1));
        img2 = load_img(hindx(2));

        hpoints = corner(img1(:,end-overlap+1:end), method, numberCorners);
        hpoints(:,1) = hpoints(:,1) + (img_size(2)-overlap);

        metric1 = cornermetric(img1, method);
        metric1 = metric1(sub2ind(img_size, hpoints(:,2), hpoints(:,1)));

        hpoints2 = corner(img2(:,1:overlap), method, numberCorners);
        metric2 = cornermetric(img2, method);
        metric2 = metric2(sub2ind(img_size, hpoints2(:,2), hpoints2(:,1)));

        all_hpoints{j} = [hpoints metric1 j*ones(size(metric1))];
        all_hpoints2{j} = [hpoints2 metric2 j*ones(size(metric2))];
      end

      tmp1 = cat(1, all_hpoints{:});
      tmp2 = cat(1, all_hpoints2{:});

      if (isempty(tmp1) || all(isnan(tmp1(:))))
        continue;
      end

      [tmp_v, tmp_indxs] = sort(tmp1(:,3), 'descend');
      tmp1 = tmp1(tmp_indxs(1:min(numberCorners, end)),:);

      [hshiftx, hshifty] = ShiftByCornerClustering(tmp1(:,1:2), tmp2(:,1:2), method, numberCorners, 1);

      for j=1:length(rindx)
        curr_pts = (tmp1(:,4) == j);

        if (any(curr_pts))

          hindx = [curr_indx(rindx(j), cindx(j)) curr_indx(rindx(j), cindx(j)+1)];

          img1 = load_img(hindx(1));
          img2 = load_img(hindx(2));

          hpoints = tmp1(curr_pts, 1:2);

          hpoints2 = LKTracker(img1, img2, hpoints, [hshiftx hshifty]);
          indices2 = CheckPointsAndNAN(img2, hpoints2);
          hpoints2 = hpoints2(indices2,:);
          hpoints = hpoints(indices2,:);

          all_hpoints{j} = hpoints;
          all_hpoints2{j} = hpoints2;
        else
          all_hpoints{j} = NaN(0,2);
          all_hpoints2{j} = NaN(0,2);
        end
      end

      hpoints = cat(1, all_hpoints{:});
      hpoints2 = cat(1, all_hpoints2{:});

      if (isempty(hpoints) || all(isnan(hpoints(:))))
        continue;
      end

      hshift = mean(hpoints - hpoints2, 1);

      if (abs(hshift(2)/hshift(1)) > prop_thresh)
        continue;
      end

      all_shifts(i,3:4) = hshift;
      all_shifts(i,6) = corr2;
      all_shifts(i,8) = size(hpoints,1);
    else
      all_shifts(i,6) = 0;
      all_shifts(i,8) = 0;
    end
  end

  % Could also use the sum/average of correlations
  [val, best] = max(sum(all_shifts(:,7:8), 2));

  if (~isnan(val))
    all_shifts(isnan(all_shifts)) = 0;
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

    img = double(imnorm(img));
  end
end

function [corrs] = correl_images(img1, img2, strength1, strength2, thresh, wthresh)

  maxs = [];

  img1 = bsxfun(@minus, img1, mean(img1, 1));
  img2 = bsxfun(@minus, img2, mean(img2, 1));

  sum1 = sum(img1.^2, 1);
  sum2 = sum(img2.^2, 1);

  nil1 = (sum1 == 0);
  nil2 = (sum2 == 0);

  strength1(nil1) = 0;
  strength2(nil2) = 0;

  goods1 = (strength1 > wthresh);
  goods2 = (strength2 > wthresh);

  strength1 = strength1 / max(strength1);
  strength2 = strength2 / max(strength2);

  indx1 = find(goods1(1:end), 1, 'last');
  indx2 = find(goods2(1:end), 1, 'first');

  if (isempty(indx1) || isempty(indx2))
    corrs = zeros(size(strength1));
    return;
  end

  corr1 = strength2.*sum(bsxfun(@times, img1(:,indx1), img2), 1) ./ sqrt(sum1(indx1) * sum2);
  corr2 = strength1.*sum(bsxfun(@times, img1, img2(:,indx2)), 1) ./ sqrt(sum1 * sum2(indx2));

  corr1 = [corr1(indx1+1:end) corr1(1:indx1)];
  corr2 = [corr2(indx2:end) corr2(1:indx2-1)];

  corr1(nil2) = 0;
  corr2(nil1) = 0;

  corrs = corr1 + corr2(end:-1:1);

  return;
end
