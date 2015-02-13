function cells = linear_movement(cells, curr_time, opts)

  dm = opts.cell_speed * opts.dt;

  cells(:, 1:2) = cells(:, 1:2) + bsxfun(@plus, opts.movement_params * dm, randn(size(cells, 1), 2) * opts.cell_variation);

  return;
end
