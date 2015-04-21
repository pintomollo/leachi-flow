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
    interp.speeds = speeds;
    interp.parameters = [vects cross lens];

    if (isstruct(opts.cell_speed))
      interp.scaling = opts.cell_speed.scaling;
    else
      interp.scaling = opts.cell_speed;
    end

    opts.cell_speed = interp;
  end

  speeds = opts.cell_speed.speeds * opts.cell_speed.scaling;
  params = opts.cell_speed.parameters.';
  origin = opts.cell_speed.segments(:,1:2).';

  dists = bsxfun(@times, (bsxfun(@plus, ...
                            bsxfun(@times, cells(:,1), params(2, :)) - ...
                              bsxfun(@times, cells(:,2), params(1,:)), ...
                            params(3,:))).^2, ...
                         params(4,:));

  frac = bsxfun(@times, bsxfun(@minus, cells(:,1), ...
                                       origin(1,:)), ...
                        params(1,:) .* params(4,:)) + ...
         bsxfun(@times, bsxfun(@minus, cells(:,2), ...
                                       origin(2,:)), ...
                        params(2,:) .* params(4,:));

  dists(frac < 0 | frac > 1) = Inf;

  [dist, indx] = min(dists, [], 2);

  %{
  figure;hold on;
  colors = jet(size(dists, 2));
  for i=1:size(dists, 2)
    goods = (indx==i);
    scatter(cells(goods, 1), cells(goods, 2), [], colors(i,:));
    plot(opts.cell_speed.segments(i, [1 3]), opts.cell_speed.segments(i, [2 4]), 'Color', colors(i,:));
  end

  keyboard
  %}

  dm = opts.cell_speed.scaling  * opts.dt * sin(2*pi*curr_time/opts.movement_params(1));
  dvect = bsxfun(@times, opts.cell_speed.speeds(indx,:) * dm, 1 + randn(size(cells, 1), opts.ndims) * opts.cell_variation);

  return;
end
