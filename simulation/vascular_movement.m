function [dvect, dt] = vascular_movement(cells, curr_time, opts)

  if (nargout == 1)
    centers = opts.creation_params.center;
    x1 = centers(1:3:end, 1);
    x2 = centers(2:3:end, 1);
    y1 = centers(1:3:end, 2);
    y2 = centers(2:3:end, 2);

    vects = [x2-x1 y2-y1];
    lens = 1 ./ sum(vects.^2, 2);
    cross = (x2.*y1 - y2.*x1);

    speeds = bsxfun(@times, vects, sqrt(lens));

    interp = get_struct('interpolation');
    interp.segments = [x1 y1 x2 y2];
    interp.speeds = [speeds; opts.creation_params.junction.vector.'];
    interp.dists = 1 ./ (opts.creation_params.property(:)).^2;
    interp.parameters = [vects cross lens sqrt(lens .* interp.dists)];

    if (isstruct(opts.cell_speed))
      interp.scaling = opts.cell_speed.scaling;
    else
      interp.scaling = opts.cell_speed;
    end

    joint = get_struct('junction');
    joint.center = opts.creation_params.junction.polygon(1:5:end,:).';
    joint.threshold = opts.creation_params.junction.threshold;
    polygon = opts.creation_params.junction.polygon;
    polygon(5:5:end,:) = [];
    polygon = cat(3, reshape(polygon(:,1), 4, []).', reshape(polygon(:,2), 4, []).');
    joint.polygon = polygon(:,[1:end 1], :);
    joint.vector = [length(x1) size(polygon, 1) max(opts.creation_params.property)^2];

    interp.joints = joint;

    opts.cell_speed = interp;

    dvect = opts;
    return;
  end

  dm = sin(2*pi*curr_time/opts.movement_params(1));
  speeds = opts.cell_speed.speeds * opts.cell_speed.scaling * dm;
  params = opts.cell_speed.parameters.';
  origin = opts.cell_speed.segments(:,1:2).';
  juncts = opts.cell_speed.joints;
  widths = opts.cell_speed.dists.';

  dists = bsxfun(@times, (bsxfun(@plus, ...
                            bsxfun(@times, cells(:,1), params(2, :)) - ...
                              bsxfun(@times, cells(:,2), params(1,:)), ...
                            params(3,:))).^2, ...
                         params(4,:) .* widths);

  frac = bsxfun(@times, bsxfun(@minus, cells(:,1), ...
                                       origin(1,:)), ...
                        params(1,:) .* params(4,:)) + ...
         bsxfun(@times, bsxfun(@minus, cells(:,2), ...
                                       origin(2,:)), ...
                        params(2,:) .* params(4,:));

  perps = bsxfun(@times, bsxfun(@minus, cells(:,1), ...
                                       origin(1,:)), ...
                        -params(2,:) .* params(5,:)) + ...
         bsxfun(@times, bsxfun(@minus, cells(:,2), ...
                                       origin(2,:)), ...
                        params(1,:) .* params(5,:));


  %dists(frac < 0 | frac > 1) = 1 + dists(frac < 0 | frac > 1);
  dists(frac < 0 | frac > 1) = Inf;

  [dist, indx] = min(dists, [], 2);

  perp = perps(sub2ind(size(perps),[1:size(perps,1)].', indx));

  indx = indx .* sign(perp);

  centering = [-params(2, abs(indx)).*sqrt(params(4, abs(indx))); params(1, abs(indx)).*sqrt(params(4, abs(indx)))];
  perp(abs(perp)<1) = 0;
  perp = perp * 0.5;

  centering = bsxfun(@times, centering, perp.');

  juncs = bsxfun(@le, bsxfun(@minus, cells(:,1), juncts.center(1, :)).^2 + ...
                        bsxfun(@minus, cells(:,2), juncts.center(2, :)).^2, ...
                      juncts.threshold);

  for i=1:juncts.vector(2)
    goods = juncs(:,i);
    if (any(goods))
      subs = cells(goods,:);
      ins = inpolygon(subs(:,1), subs(:,2), juncts.polygon(i,:,1), juncts.polygon(i,:,2));
      goods(goods) = ins;

      indx(goods) = i+juncts.vector(1);
    end
  end

  %dvect = bsxfun(@times, (opts.cell_speed.speeds(abs(indx),:) + centering.') * dm, 1 + randn(size(cells, 1), opts.ndims) * opts.cell_variation);
  dvect = bsxfun(@times, (speeds(abs(indx),:) - centering.'), 1 + randn(size(cells, 1), opts.ndims) * opts.cell_variation);

  %{
  figure;hold on;
  ncols = size(opts.cell_speed.speeds, 1);
  colors = jet(4*ncols);
  for i=1:ncols
    goodsp1 = (indx==i & perp > 1);
    goodsp2 = (indx==i & perp <= 1);
    goodsn1 = (indx==-i & perp < -1);
    goodsn2 = (indx==-i & perp >= -1);
    scatter(cells(goodsp1, 1), cells(goodsp1, 2), [], colors(4*i-3,:));
    scatter(cells(goodsp2, 1), cells(goodsp2, 2), [], colors(4*i-2,:));
    scatter(cells(goodsn1, 1), cells(goodsn1, 2), [], colors(4*i-1,:));
    scatter(cells(goodsn2, 1), cells(goodsn2, 2), [], colors(4*i,:));

    if (i > juncts.vector(1))
      plot(juncts.polygon(i-juncts.vector(1),:, 1), juncts.polygon(i-juncts.vector(1),:, 2), 'Color', colors(4*i,:));
    else
      plot(opts.cell_speed.segments(i, [1 3]), opts.cell_speed.segments(i, [2 4]), 'Color', colors(4*i,:));
    end
  end
  quiver(cells(:,1), cells(:,2), dvect(:,1), dvect(:,2))
  plot(opts.creation_params.border(:,1), opts.creation_params.border(:,2), 'k');
  keyboard
  %}

  dlen = sum(dvect.^2, 2);
  dt = min(opts.dmax/max(dlen), opts.dt);

  dvect = dvect * dt;

  return;
end
