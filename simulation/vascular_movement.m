function [dvect, opts] = vascular_movement(cells, curr_time, opts)

  if (nargout > 1)
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
    interp.parameters = [vects cross lens];
    interp.dists = 1 ./ (opts.creation_params.property(:).').^2;

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
  end

  speeds = opts.cell_speed.speeds * opts.cell_speed.scaling;
  params = opts.cell_speed.parameters.';
  origin = opts.cell_speed.segments(:,1:2).';
  juncts = opts.cell_speed.joints;
  widths = opts.cell_speed.dists;

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

  perp = bsxfun(@times, bsxfun(@minus, cells(:,1), ...
                                       -origin(2,:)), ...
                        params(1,:) .* widths) + ...
         bsxfun(@times, bsxfun(@minus, cells(:,2), ...
                                       origin(1,:)), ...
                        params(2,:) .* widths);


  dists(frac < 0 | frac > 1) = 1 + dists(frac < 0 | frac > 1);

  [dist, indx] = min(dists, [], 2);

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

  %{
  figure;hold on;
  ncols = size(opts.cell_speed.speeds, 1);
  colors = jet(ncols);
  for i=1:ncols
    goods = (indx==i);
    scatter(cells(goods, 1), cells(goods, 2), [], colors(i,:));

    if (i > juncts.vector(1))
      plot(juncts.polygon(i-juncts.vector(1),:, 1), juncts.polygon(i-juncts.vector(1),:, 2), 'Color', colors(i,:));
    else
      plot(opts.cell_speed.segments(i, [1 3]), opts.cell_speed.segments(i, [2 4]), 'Color', colors(i,:));
    end
  end
  keyboard
  %}

  dm = opts.cell_speed.scaling  * opts.dt * sin(2*pi*curr_time/opts.movement_params(1));
  dvect = bsxfun(@times, opts.cell_speed.speeds(indx,:) * dm, 1 + randn(size(cells, 1), opts.ndims) * opts.cell_variation);

  return;
end
