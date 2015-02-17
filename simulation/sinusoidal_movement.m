function cells = quadratic_movement(cells, curr_time, opts)

  dm = opts.cell_speed * opts.dt;
  ncells = size(cells, 1);

  if (size(opts.movement_params, 1) > 1)
    intens = [[sin((cells(:,1) - opts.movement_params(1,1)) / opts.movement_params(1,2)) + ...
               sin((cells(:,2) - opts.movement_params(1,3)) / opts.movement_params(1,4))], ...
              [sin((cells(:,1) - opts.movement_params(2,1)) / opts.movement_params(2,2)) + ...
               sin((cells(:,2) - opts.movement_params(2,3)) / opts.movement_params(2,4))]];
  else
    intens = [sin((cells(:,1) - opts.movement_params(1)) / opts.movement_params(2)) + ...
              sin((cells(:,2) - opts.movement_params(3)) / opts.movement_params(4))];
  end

  cells(:, 1:2) = cells(:, 1:2) + bsxfun(@times, intens * dm, 1 + (randn(ncells, 2)*opts.cell_variation));

  %cells(:, 1:2) = cells(:, 1:2) + bsxfun(@plus, intens * dm, randn(size(cells, 1), 2) * opts.cell_variation);
  %cells(:, 1:2) = cells(:, 1:2) + bsxfun(@plus, intens * dm, ones(ncells, 2));

  return;
end
