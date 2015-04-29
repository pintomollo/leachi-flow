function cells = uniform_init(opts)

  density = opts.cell_density * opts.pixel_size^2;

  curr_size = (opts.image_size-1) * (1 + 2*opts.outside_ridge);
  curr_n = ceil(density * prod(curr_size));

  b_leachi = get_struct('botrylloides_leachi');
  cell_props = gmdistribution(b_leachi.blood_cell.mu, b_leachi.blood_cell.sigma, b_leachi.blood_cell.proportions);

  pos = bsxfun(@minus, bsxfun(@times, rand(curr_n, opts.ndims), curr_size) + 1, ...
                       (opts.image_size-1)*opts.outside_ridge);

  cells = [pos random(cell_props, curr_n)];
  cells = cells(~any(cells(:,3:4)<=0, 2), :);

  return;
end
