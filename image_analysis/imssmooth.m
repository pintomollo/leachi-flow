function img = imssmooth(img, sigma, thresh)

  if (nargin == 1)
    sigma = 0.5;
    thresh = 0;
  elseif (nargin == 2)
    thresh = 0;
  end

  ssizes = size(img);
  ndims = length(ssizes);

  if (isa(img, 'ndSparse'))
    dims = [2:ndims 1];
    for i=1:ndims
      img = sparse(img);
      img = gaussian_sparse_mex(img, sigma, thresh);
      img = ndSparse(img, ssizes);
      img = permute(img, dims);
      ssizes = ssizes(dims);
    end
  else
    img = gaussian_sparse_mex(img, sigma, thresh);
    img = permute(img, [2 1]);
    img = gaussian_sparse_mex(img, sigma, thresh);
    img = permute(img, [2 1]);
  end

  return;
end
