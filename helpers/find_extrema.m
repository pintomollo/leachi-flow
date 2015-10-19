function [xmax, imax] = find_extrema(x, nneigh)

  if (nargin == 1 || isempty(nneigh))
    nneigh = 5;
  end

  ssize = size(x);
  ndim = length(ssize);

  if (ndim == 2 && (ssize(1) == 1 || ssize(2) == 1))
    x = x(:);
    ndim = 1;
  end

  x(isnan(x)) = min(x(:));

  nneigh = [nneigh(:).' ones(1, ndim-length(nneigh))*nneigh(1)];
  nneigh = nneigh(1:ndim);

  if (ndim == 1)
    [xmax, imax] = local_extrema(x, nneigh);
  else
    [xmax, imax] = nd_extrema(x, nneigh);

    cols = cell(ndim, 1);
    [cols{:}] = ind2sub(ssize, imax);
    imax = cat(2, cols{:});
  end

  return;
end

function [xmax, imax] = nd_extrema(x, nneigh)

  minval = min(x(:));

  nneigh = 2*nneigh;
  kernel = true(nneigh+1);

  maxval = imdilate(x, kernel);
  maxpos = (x == maxval & x ~= minval);

  maxval = maxval(maxpos);

  if (isempty(maxval))
    xmax = minval;
    imax = ones(1, ndims(x));
  end

  bw = imdilate(maxpos, kernel);

  labels = labelmatrix(bwconncomp(bw));
  labels = double(labels(maxpos));

  maxdis = bwdist(~bw);
  maxdis(maxpos) = maxdis(maxpos) + rand(sum(maxpos(:)), 1);
  maxdis = maxdis(maxpos);

  maxpos = find(maxpos);

  [maxval, indxs] = sortrows([labels -maxval -maxdis]);
  firsts = (diff([0; maxval(:, 1)]) > 0);

  imax = maxpos(indxs(firsts));
  xmax = -maxval(firsts, 2);

  return;
end
