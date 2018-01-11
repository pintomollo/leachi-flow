function [hfig, nelems] = density_scatter(varargin)

  haxes = [];
  x = [];
  y = [];
  resol = 1;
  map = @(x)(brewermap(x, '*RdGy'));

  if (ishandle(varargin{1}))
    haxes = varargin{1};
    varargin(1) = [];
  end

  nargin = length(varargin);
  if (nargin == 1)
    x = varargin{1};

    y = x(:,2);
    x = x(:,1);
  elseif (nargin ==2)
    x = varargin{1};
    y = varargin{2};

    if (isa(y, 'function_handle'))
      map = y;
      y = x(:,2);
      x = x(:,1);
    elseif (numel(y) ~= numel(x))
      resol = y;
      y = x(:,2);
      x = x(:,1);
    end
  elseif (nargin ==3)
    x = varargin{1};
    y = varargin{2};
    resol = varargin{3};

    if (isa(resol, 'function_handle'))
      map = resol;
    end
  elseif (nargin ==4)
    x = varargin{1};
    y = varargin{2};
    resol = varargin{3};
    map = varargin{4};
  end

  if (numel(resol) < 2)
    resol = resol([1 1]);
  end

  nbins = ceil([range(x) max(abs(y))] ./ resol);

  xbins = [0:nbins(1)];
  ybins = [-nbins(2)-0.5:nbins(2)+0.5];

  xbins = xbins * (range(x)/xbins(end));
  ybins = ybins * (max(abs(y))/ybins(end));

  xbins = xbins + min(x);

  xc = diff(xbins) + xbins(1:end-1);
  yc = diff(ybins) + ybins(1:end-1);

  xbins(1) = xbins(1)-1;
  xbins(end) = xbins(end)+1;
  ybins(1) = ybins(1)-1;
  ybins(end) = ybins(end)+1;

  [density, pos_edges, xs, pos] = histcn([x y], xbins, ybins);
  [xres, yres, vals] = find(sparse(density));
  [v, junk, indxj] = unique(vals);
  [junk, indxr] = sort(indxj);

  nv = max(v);
  colors = map(nv);
  colv = colors(v,:);
  colv = colv(indxj,:);
  colv = colv(indxr,:);

  xres = xc(xres);
  yres = yc(yres);

  xres = xres(indxr);
  yres = yres(indxr);
  vres = vals(indxr);

  if (~isempty(haxes))
    h = scatter(haxes, xres, yres, [], colv, 'filled');
  else
    h = scatter(xres, yres, [], colv, 'filled');
    haxes = get(h, 'Parent');
  end
  set(get(h, 'Child'), 'CDataMapping', 'direct');
  colormap(haxes, colors);
  caxis([1 nv]);

  if (nargout > 0)
    hfig = h;

    if (nargout > 1)
      nelems = length(xres);
    end
  end

  return;
end
