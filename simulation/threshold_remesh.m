function [cells, thresh] = threshold_remesh(cells, dvect, thresh, opts)

  % The dimensionaltiy of our problem
  dims = [1:opts.ndims];

  % The next position
  cells(:, dims) = cells(:, dims) + dvect;

  % The fraction of the threshold traveled by the fastest particle
  dthresh = sqrt(max(sum(dvect.^2, 2))) / (opts.outside_ridge * max(opts.image_size));

  % If we have traveled all the threshold, we need ot update the particles
  thresh = thresh + dthresh;
  if (thresh > opts.remesh_params)

    % Count the cells outside of the image but still inside the threshold
    rim = (opts.image_size-1) * opts.outside_ridge * (1 - opts.remesh_params);
    rimmed =  all(bsxfun(@gt, cells(:,dims), -rim) & bsxfun(@lt, cells(:,dims), opts.image_size+rim), 2) & ...
             ~all(bsxfun(@gt, cells(:,dims), [1 1]) & bsxfun(@lt, cells(:,dims), opts.image_size), 2);

    % Check which ones enter and wich ones leave the image
    center = opts.image_size/2;
    pos = sum(bsxfun(@minus, cells(rimmed, dims), center).^2, 2);
    prev = sum(bsxfun(@minus, cells(rimmed, dims) - dvect(rimmed, :), center).^2, 2);
    inners = (pos < prev);

    % Update the density of create particles accordingly to get a more or less constant density inside the image
    if (any(inners) && ~all(inners))
      ratio = sum(inners);
      ratio = ratio / (length(pos) - ratio);
      opts.cell_density = opts.cell_density / (0.2*(3 + 2*ratio));
    end

    % Create the new cells (throughout the image actually)
    new_cells = opts.create_cells(opts);

    % Keep the current valid ones and the new outsiders
    valids = all(bsxfun(@gt, cells(:,dims), -rim) & bsxfun(@lt, cells(:,dims), opts.image_size+rim), 2);
    new_valids = all(bsxfun(@gt, new_cells(:,dims), -rim) & bsxfun(@lt, new_cells(:,dims), opts.image_size+rim), 2);

    % Merge them together
    cells = [new_cells(~new_valids, :); cells(valids, :)];

    thresh = 0;
  end

  return;
end
