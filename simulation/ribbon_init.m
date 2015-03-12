function cells = ribbon_init(opts)

  img_size = opts.image_size - 1;
  center = img_size / 2;
  angle = opts.creation_params(1);

  density = opts.n_cells / (min(img_size)*opts.creation_params(2));
  curr_size = [sqrt(sum(img_size.^2)) * (1 + 2*opts.outside_ridge) opts.creation_params(2)];
  curr_n = ceil(density * prod(curr_size));

  pos = bsxfun(@times, rand(curr_n, opts.ndims)-0.5, curr_size);
  rot = [pos(:,1)*cos(angle) - pos(:,2)*sin(angle), ...
         pos(:,1)*sin(angle) + pos(:,2)*cos(angle)];

  pos = bsxfun(@plus, rot, center);
  cells = [pos randn(curr_n, opts.nprops)*opts.cell_variation + opts.cell_size];

  return;
end
