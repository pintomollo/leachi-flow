function [img, background] = draw_background(img_size, background, curr_time, opts)

  img = zeros(img_size);
  if (~isempty(background))
    if (size(background.vessel, 2) == 2)
      x1 = background.vessel(1:end-1, 1);
      x2 = background.vessel(2:end, 1);
      y1 = background.vessel(1:end-1, 2);
      y2 = background.vessel(2:end, 2);

      valids = ~(isnan(x1 + y1 + x2 + y2));
      x1 = x1(valids);
      x2 = x2(valids);
      y1 = y1(valids);
      y2 = y2(valids);

      pts = unique([x1 y1; x2 y2], 'rows');

      vects = [x2-x1 y2-y1];
      lens = 1 ./ sum(vects.^2, 2);
      cross = (x2.*y1 - y2.*x1);
      origin = [x1 y1].';

      params = [vects cross lens sqrt(lens)].';

      [X,Y] = meshgrid([1:img_size(2)], [1:img_size(1)]);

      pixels = [X(:) Y(:)];

      dists = bsxfun(@times, (bsxfun(@plus, ...
                                bsxfun(@times, pixels(:,1), params(2, :)) - ...
                                  bsxfun(@times, pixels(:,2), params(1,:)), ...
                                params(3,:))).^2, ...
                             params(4,:));

      knots = bsxfun(@minus, pixels(:,1), pts(:,1).').^2 + ...
              bsxfun(@minus, pixels(:,2), pts(:,2).').^2;

      frac = bsxfun(@times, bsxfun(@minus, pixels(:,1), ...
                                           origin(1,:)), ...
                            params(1,:) .* params(4,:)) + ...
             bsxfun(@times, bsxfun(@minus, pixels(:,2), ...
                                           origin(2,:)), ...
                            params(2,:) .* params(4,:));

      knots = repmat(min(knots, [], 2), 1, size(frac,2));

      dists(frac > 1 | frac < 0) = knots(frac > 1 | frac < 0);
      dists = min(dists, [], 2);

      img = reshape(exp(-dists/(2*background.vessel_properties(2)^2))*background.vessel_properties(1), img_size);

      if (~isempty(background.clot))
        img = img + draw_gaussians_mex(img_size, background.clot);
      end

      background.vessel = img;
    else
      img = background.vessel;
    end

    if (~isempty(background.ampulla) && nargin > 2)
      variation = 0.2;
      dm = 1 + variation * cos(2*pi*curr_time/opts.movement_params(1));
      curr_size = background.ampulla(:,3) * dm;
      img = img + draw_bumps_mex(img_size, [background.ampulla(:, 1:2) curr_size background.ampulla(:, 4)], 0.2);
    end
  end

  return;
end
