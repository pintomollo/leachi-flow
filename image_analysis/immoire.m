function [img, moire] = immoire(img, thresh)

  if (nargin < 2)
    thresh = 5;
  end

  F = fft2(img);

  outliers1 = filter_thresh(real(F), thresh);
  outliers2 = filter_thresh(imag(F), thresh);

  F1 = F;
  F2 = F;
  F3 = F;

  F1(outliers1) = 0;
  F2(outliers2) = 0;
  F3(outliers1 | outliers2) = 0;

  G1 = real(ifft2(F1));
  G2 = real(ifft2(F2));
  G3 = real(ifft2(F3));

  figure;
  subplot(2,3,1);imagesc(G1);
  subplot(2,3,2);imagesc(G2);
  subplot(2,3,3);imagesc(G3);

  subplot(2,3,4);imagesc(G1-img);
  subplot(2,3,5);imagesc(G2-img);
  subplot(2,3,6);imagesc(G3-img);

  keyboard

  return;
end

function res = filter_thresh(data, thresh)

  res = false(size(data));
  data = data(:);

  % We need to robustly estimate the mean and the standard deviation of the pixel
  % distribution. In opposition with what is proposed in [1], we utilize the median
  % and the MAD estimators for this purpose.
  dmed = median(data);
  dmad = 1.4826 * median(abs(data-dmed));

  % Here we most likely have a problem, so try another approach to estimate the
  % standard deviation.
  if (dmad==0)
    dmad = 1.4826 * mean(abs(data-dmed));

    % If nothing else worked, follow [1] to esimate the standard deviation
    if (dmad == 0)
      dmad = sqrt(mean(data.^2 - dmed^2));
    end
  end

  % Sort the data, keeping track of the correspondency
  [data, indx] = sort(data-dmed);

  figure;hist(data, 200);

  % And measure the distance between pixel values
  dist = [0; diff(data)];

  % Our cosmic rays are higher than the estimated mean and have a gap larger than
  % the standard deviation times the threshold
  bads = (data > 0 & dist > thresh*dmad);

  % If we find one, [1] removes all the higher-valued pixels. Given that we have
  % ordere them, we find the first index and replace all higher values by the median
  if (any(bads))
    first = find(bads, 1, 'first');
    res(indx(first:end)) = true;
  end

  % Our cosmic rays are higher than the estimated mean and have a gap larger than
  % the standard deviation times the threshold
  bads = (data < 0 & dist > thresh*dmad);

  % If we find one, [1] removes all the higher-valued pixels. Given that we have
  % ordere them, we find the first index and replace all higher values by the median
  if (any(bads))
    first = find(bads, 1, 'last');
    res(indx(1:first)) = true;
  end

  return;
end
