function cells = uniform_init(simul, opts)

  density = simul.cell_density * opts.pixel_size^2;

  curr_size = (simul.image_size-1) * (1 + 2*simul.outside_ridge);
  curr_n = ceil(density * prod(curr_size));

  b_leachi = get_struct('botrylloides_leachi');
  cell_props = gmdistribution(b_leachi.blood_cell.mu, b_leachi.blood_cell.sigma, b_leachi.blood_cell.proportions);

  pos = bsxfun(@minus, bsxfun(@times, rand(curr_n, simul.ndims), curr_size) + 1, ...
                       (simul.image_size-1)*simul.outside_ridge);

  cells = [pos random(cell_props, curr_n)];
  cells = cells(~any(cells(:,3:4)<=0, 2), :);

  return;
end
