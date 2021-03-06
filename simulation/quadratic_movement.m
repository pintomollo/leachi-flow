function dvect = quadratic_movement(cells, curr_time, opts)

  dm = opts.cell_speed * opts.dt;
  ncells = size(cells, 1);

  grid = [cells(:,1).^2, cells(:,1).*cells(:,2), cells(:,2).^2 ones(ncells, 1)];
  if (size(opts.movement_params, 1) > 1)
    intens = [sum(bsxfun(@times, grid, opts.movement_params(1, 1:4)), 2), ...
              sum(bsxfun(@times, grid, opts.movement_params(2, 1:4)), 2)];
  else
    intens = sum(bsxfun(@times, grid, opts.movement_params(1:4)), 2);
  end

  intens = imnorm(abs(intens)).*sign(intens);

  dvect = bsxfun(@times, intens * dm, 1 + (randn(ncells, opts.ndims)*opts.cell_variation));

  return;
end
