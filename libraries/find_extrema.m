function find_extrema(x, nneigh)

  if (nargin == 1 || isempty(nneigh))
    nneigh = 5;
  end

  ssize = size(x);
  ndim = length(ssize);

  if (ndim == 2 && (ssize(1) == 1 || ssize(2) == 1))
    x = x(:);
    ndim = 1;
  end

  nneigh = [nneigh(:).' ones(1, ndim-length(nneigh))*nneigh(1)];
  nneigh = nneigh(1:ndim);

  if (ndim == 1)
    kernel = true(2*nneigh+1, 1);
    shrink = true(3, 1);
  else
    kernel = true(2*nneigh+1);
    shrink = true(ones(1, ndim)*3);
  end

  maxval = imdilate(x, kernel);
  maxpos = (x == maxval);

  maxunique = myshrink(maxpos, kernel, shrink);

  keyboard

  return;
end

function pts = myshrink(bw, dil, ero)

  %{  %% From algbwmorph
      case 'shrink'
        
        lut   = [];
        table = images.internal.lutshrink();
        
        % First subiteration
        m   = bwlookup(bw, table);
        sub = bw & ~m;
        bw(1:2:end,1:2:end) = sub(1:2:end,1:2:end);
        
        % Second subiteration
        m   = bwlookup(bw, table);
        sub = bw & ~m;
        bw(2:2:end,2:2:end) = sub(2:2:end,2:2:end);
        
        % Third subiteration
        m   = bwlookup(bw, table);
        sub = bw & ~m;
        bw(1:2:end,2:2:end) = sub(1:2:end,2:2:end);
        
        % Fourth subiteration
        m   = bwlookup(bw, table);
        sub = bw & ~m;
        bw(2:2:end,1:2:end) = sub(2:2:end,1:2:end);
  %}

  pts = false(size(bw));
  niter = sum(bw(:));

  for i=1:niter
    step = imerode(bw, ero);
    areas = imdilate(step, dil);

    pts = (pts | (bw & (~areas)));

    if (~any(step(:)) || ~any(xor(step(:),bw(:))))
      break;
    end

    bw = step;
  end

  return;
end
