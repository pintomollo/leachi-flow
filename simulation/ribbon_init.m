function cells = ribbon_init(opts)

  b_leachi = get_struct('botrylloides_leachi');
  cell_props = gmdistribution(b_leachi.blood_cell.mu, b_leachi.blood_cell.sigma, b_leachi.blood_cell.proportions);

  img_size = opts.image_size - 1;
  center = img_size / 2;
  angle = opts.creation_params(1);

  density = opts.cell_density * opts.pixel_size^2;
  curr_size = [sqrt(sum(img_size.^2)) * (1 + 2*opts.outside_ridge) opts.creation_params(2)];
  curr_n = ceil(density * prod(curr_size));

  pos = bsxfun(@times, rand(curr_n, opts.ndims)-0.5, curr_size);
  rot = [pos(:,1)*cos(angle) - pos(:,2)*sin(angle), ...
         pos(:,1)*sin(angle) + pos(:,2)*cos(angle)];

  pos = bsxfun(@plus, rot, center);

  cells = [pos random(cell_props, curr_n)];
  cells = cells(~any(cells(:,3:4)<=0, 2), :);
  cells(:,3) = cells(:,3) / opts.pixel_size;

  return;
end
