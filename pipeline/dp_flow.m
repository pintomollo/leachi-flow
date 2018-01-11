function [x, y] = dp_flow(t, v, resol)

  if (nargin == 1)
    v = t(:,2);
    t = t(:,1);
    resol = 1;
  elseif (nargin ==2)
    if (numel(v)==numel(t))
      resol = 1;
    else
      resol = v;
      v = t(:,2);
      t = t(:,1);
    end
  end

  if (numel(resol) < 2)
    resol = resol([1 1]);
  end

  %nbins = nbins + (1 - mod(nbins, 2));
  nbins = ceil([range(t) max(abs(v))] ./ resol);

  xbins = [0:nbins(1)];
  ybins = [-nbins(2)-0.5:nbins(2)+0.5];

  xbins = xbins * (range(t)/xbins(end));
  ybins = ybins * (max(abs(v))/ybins(end));

  xbins = xbins + min(t);

  %disp([xbins(2)-xbins(1) ybins(2)-ybins(1)])

  xbins(1) = xbins(1)-1;
  xbins(end) = xbins(end)+1;
  ybins(1) = ybins(1)-1;
  ybins(end) = ybins(end)+1;

  [map, pos_edges, xs, pos] = histcn([t v], xbins, ybins);

  p = get_struct('parameter');
  p.weights.alpha = 0.45;
  p.weights.beta = 0.9;
  p.weights.gamma = 1/9;
  weight_flow(map, p.weights);

  p.parameters.nhood   = 50;      % Neighborhood size
  p.parameters.force_circularity = false;
  p.parameters.dp_method = 'normal';
  p.parameters.alpha   = 0.35; % Prop. of smoothness VS data
  p.parameters.beta    = 0.85; % Prop. of path VS intensity
  p.parameters.gamma   = 0.25; % Prop. of dx VS d2x

  nhood = max(ceil(length(ybins)/5), 7);
  p.parameters.init = 1-imnorm(differentiator(cumsum(sum(map(1:5,:))), 1, nhood));
  p.parameters.final = find(imnorm(differentiator(cumsum(sum(map(1:5,:))), 1, nhood))>0.75);

  orig_paths = dynamic_programming(1-imnorm(map), p.parameters, @weight_flow, p.weights, []);

  x = diff(xbins) + xbins(1:end-1);
  y = diff(ybins) + ybins(1:end-1);

  y = y(orig_paths);

  %{
  %vars = [0:0.05:1];
  vars = [1:20];
  colors=redbluemap(length(vars));
  for i=1:length(vars)
    %p.weights.alpha = vars(i);
    p.weights.gamma = 1/vars(i);

    paths = dynamic_programming(1-imnorm(map), p.parameters, @weight_flow, p.weights, []);

    r = y(paths);

    plot(x, r, 'Color', colors(i,:))
  end
  plot(x, y(orig_paths), 'g')
  %}

  %[pp] = csaps(x, y, 1/(x(end)));

  %figure;imagesc(map)
  %hold on;
  %plot(paths, 1:length(paths), 'k');

  %figure;hold on;
  %plot(x, y);
  %fnplt(pp, 'k')

  %keyboard

  return;
end

function w = weight_flow(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = 1/params.gamma;

  img = imnorm(img);
  [nrows,npts] = size(img);

  center = gaussian_mex(img, gamma);

  half = (npts-1)/2 + 1;
  folded = 2 - (center(:,half:-1:1) + center(:,half:end));
  domain = cumsum(folded, 2);
  stats = imnorm(sum(folded, 1));
  domain = bsxfun(@minus, 2*domain, domain(:,end));

  center = imnorm(center, [], [], 'rows');
  %domain = imnorm(1-domain, [], [], 'rows');
  %domain = imnorm(1-domain, [], [], 'rows');
  domain = gaussian_mex(double(domain < 0.99), gamma);

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

  domain_value = domain_value(:,[end:-1:1 2:end]);

  w = alpha * center + (1-alpha) * domain_value;

  %dist = [1:half];
  %idist = max(npts - dist, 1);f

  %domain = bsxfun(@rdivide, domain, dist);
  %outside = bsxfun(@rdivide, outside, idist);

  %domain_value = (1-domain)*beta + (1-beta)*outside;

  %figure;imagesc(center)
  %figure;imagesc(domain)
  %figure;plot(stats)
  %figure;imagesc(outside)
  %figure;imagesc(domain_value)
  %figure;imagesc(w)

  return;
end
