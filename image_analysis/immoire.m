function [img, moire] = immoire(img, thresh, max_size)

  if (nargin < 2)
    thresh = 8;
  end
  if (nargin < 3)
    max_size = 5;
  end

  F = fft2(img);
  F = fftshift(F);

  G = log(abs(F));

  noise = estimate_noise(G);

  G = gaussian_mex(G, 0.67);

  intens_thresh = noise(1) + thresh*noise(2);

  % Create a mask to identify the local maxima
  tmp_size = ceil(max_size/2);
  tmp_size = tmp_size + mod(tmp_size+1, 2);
  mask = ones(tmp_size);
  mask((end-1)/2+1) = 0;

  % If need be, perpare the averaging mask
  mask_avg = ones(tmp_size);
  mask_avg = mask_avg / numel(mask_avg);

  % Performs the actual spot detection "a trous" algorithm [1]
  atrous = imatrou(G, max_size, thresh);

  % Compute the local average
  avgs = imfilter(G, mask_avg);

  % Get the local maxima
  bw = (atrous > 0) & (G >= imdilate(G, mask)) & (avgs >= intens_thresh);

  % Shrink them to single pixel values
  bw = bwmorph(bw, 'shrink', Inf);

  center = floor(size(img)/2)+1;
  bw(center(1), center(2)) = false;

  if (any(bw(:)))
    moire = imnorm(gaussian_mex(double(bw), max_size));
  else
    moire = double(bw);
  end

  H = F .* (1-moire);
  img = real(ifft2(ifftshift(H)));

  H = F .* moire;
  moire = real(ifft2(ifftshift(H)));

  return;
end
