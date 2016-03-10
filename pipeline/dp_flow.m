function dp_flow(t, v, nbins)

  if (nargin == 1)
    v = t(:,2);
    t = t(:,1);
    nbins = 64;
  elseif (nargin ==2)
    if (numel(v)==numel(t))
      nbins = 64;
    else
      nbins = v;
      v = t(:,2);
      t = t(:,1);
    end
  end

  if (numel(nbins) < 2)
    nbins = nbins([1 1]);
  end

  xbins = [0:nbins(1)];
  ybins = [-nbins(2)-0.5:nbins(2)+0.5];

  xbins = xbins * (range(t)/xbins(end));
  ybins = ybins * (max(abs(v))/ybins(end));

  xbins = xbins + min(t);

  xbins(1) = xbins(1)-eps;
  xbins(end) = xbins(end)+eps;
  ybins(1) = ybins(1)-eps;
  ybins(end) = ybins(end)+eps;

  [map, pos_edges, xs, pos] = histcn([t v], xbins, ybins);

  figure;imagesc(map)

  p = get_struct('parameter');
  p.weights.alpha = 0.5;
  p.weights.beta = 0.75;
  p.weights.gamma = 1/2;
  weight_flow(map, p.weights);

  p.parameters.nhood   = 5;      % Neighborhood size
  p.parameters.alpha   = 0.5; % Prop. of smoothness VS data
  p.parameters.beta    = 0.5; % Prop. of path VS intensity
  p.parameters.gamma   = 0.5; % Prop. of dx VS d2x

  paths = dynamic_programming(map, p.parameters, @weight_flow, p.weights, []);
  keyboard

  return;
end

function w = weight_flow(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = 1/params.gamma;

  img = imnorm(img);
  [nrows,npts] = size(img);

  center = gaussian_mex(1-img, gamma);

  half = (npts-1)/2 + 1;
  folded = center(:,half:-1:1) + center(:,half:end);
  domain = cumsum(folded, 2);
  stats = imnorm(sum(folded, 1));
  domain = bsxfun(@minus, 2*domain, domain(:,end));

  center = imnorm(center, [], [], 'rows');
  domain = imnorm(1-domain, [], [], 'rows');

  limit = find(stats > 0.95, 1, 'first');
  %sqrt(-2*log(0.001))
  sl0001 = 3.7169;
  %sqrt(-2*log(0.999))
  sl0999 = 0.0447;
  o = ((half-limit)/(sl0999 - sl0001));
  u = limit - o*sl0001;
  outside = exp(-([1:half] - u).^2 / (2*o^2));
  outside = repmat(outside, nrows, 1);

  domain_value = domain*beta + (1-beta)*outside;

  w = domain_value(:,[end:-1:1 2:end]);

  %dist = [1:half];
  %idist = max(npts - dist, 1);

  %domain = bsxfun(@rdivide, domain, dist);
  %outside = bsxfun(@rdivide, outside, idist);

  %domain_value = (1-domain)*beta + (1-beta)*outside;

  %figure;imagesc(center)
  %figure;imagesc(domain)
  %figure;plot(stats)
  %figure;imagesc(outside)
  %figure;imagesc(domain_value)
  figure;imagesc(w)

  return;
end
