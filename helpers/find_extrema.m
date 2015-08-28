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

  switch ndim
    case 1
      [xmax, imax] = local_extrema(x, nneigh);
    case 2
      [xmax, imax] = dd_extrema(x, nneigh);
    otherwise
      %% I was thinking to use the 
      error('not yet implemented');
      [xmax, imax] = nd_extrema(x, nneigh);
  end

  return;
end

function [xmax, imax] = dd_extrema(x, nneigh)

  nneigh = 2*nneigh;
  kernel = true(nneigh+1);

  maxval = imdilate(x, kernel);
  maxpos = (x == maxval);

  maxpos = padarray(maxpos, nneigh, false, 'both');

  maxpos = bwmorph(imclose(maxpos, kernel), 'shrink', Inf);

  maxpos = maxpos(nneigh(1)+1:end-nneigh(1), nneigh(1)+1:end-nneigh(1));

  [i,j] = find(maxpos);
  imax = [i(:), j(:)];

  xmax = x(maxpos);
  xmax = xmax(:);

  return;
end

function [xmax, imax] = nd_extrema(x, nneigh)

  kernel = true(2*nneigh+1);

  shrink = true(ones(1, ndim)*3);
  shrink(ceil(end/2)) = false;

  maxval = imdilate(x, kernel);
  bw = (x == maxval);

  pts = false(size(bw));
  niter = sum(bw(:));

  pts = imclose(pts, kernel);

  for i=1:niter
    step = imerode(bw, shrink);
    areas = imdilate(step, kernel);

    pts = (pts | (bw & (~areas)));

    if (~any(step(:)) || ~any(xor(step(:),bw(:))))
      break;
    end

    bw = step;
  end

  %% Then bwlabel, loop over them, pickup the maxs, use bwdist to chose centered pixel

  return;
end
