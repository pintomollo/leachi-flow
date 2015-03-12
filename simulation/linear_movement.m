function dvect = linear_movement(cells, curr_time, opts)

  dm = opts.cell_speed * opts.dt;
  dvect = bsxfun(@plus, opts.movement_params * dm, randn(size(cells, 1), opts.ndims) * opts.cell_variation);

  return;
end
