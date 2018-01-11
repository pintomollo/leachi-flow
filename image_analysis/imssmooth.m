function img = imssmooth(img, dim, sigma, thresh)

  if (nargin == 1)
    dim = [];
    sigma = 0.5;
    thresh = 0;
  elseif (nargin == 2)
    sigma = 0.5;
    thresh = 0;
  elseif (nargin == 3)
    thresh = 0;
  end

  ssizes = size(img);
  all_dim = [1:length(ssizes)];

  if (isempty(dim))
    ndims = length(ssizes);
    dim = [1:ndims];
  else
    ndims = length(dim);
  end

  is_nd = (isa(img, 'ndSparse'));

  for i=1:ndims
    perm_dim = all_dim;
    perm_dim(dim(i)) = 1;
    perm_dim(1) = dim(i);

    img = permute(img, perm_dim);

    if (is_nd)
      img = sparse(img);
      img = gaussian_sparse_mex(img, sigma, thresh);
      img = ndSparse(img, ssizes(perm_dim));
    else
      img = gaussian_sparse_mex(img, sigma, thresh);
    end

    img = ipermute(img, perm_dim);
  end

  return;
end
