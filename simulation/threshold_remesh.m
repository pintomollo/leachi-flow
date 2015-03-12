function [cells, thresh] = threshold_remesh(cells, dvect, thresh, opts)

  dims = [1:opts.ndims];
  cells(:, dims) = cells(:, dims) + dvect;

  dthresh = max(sum(dvect.^2, 2)) / (opts.outside_ridge * max(opts.image_size));

  thresh = thresh + dthresh;
  if (thresh > opts.remesh_params)
    new_cells = opts.create_cells(opts);

    rim = (opts.image_size-1) * opts.outside_ridge * (1 - opts.remesh_params);

    valids = all(bsxfun(@gt, cells(:,dims), -rim) & bsxfun(@lt, cells(:,dims), opts.image_size+rim), 2);
    new_valids = all(bsxfun(@gt, new_cells(:,dims), -rim) & bsxfun(@lt, new_cells(:,dims), opts.image_size+rim), 2);

    cells = [new_cells(~new_valids, :); cells(valids, :)];

    thresh = 0;
  end

  return;
end
