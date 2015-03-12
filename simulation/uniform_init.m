function cells = uniform_init(opts)

  density = opts.n_cells / prod(opts.image_size-1);

  curr_size = (opts.image_size-1) * (1 + 2*opts.outside_ridge);
  curr_n = ceil(density * prod(curr_size));


  pos = bsxfun(@minus, bsxfun(@times, rand(curr_n, opts.ndims), curr_size) + 1, ...
                       (opts.image_size-1)*opts.outside_ridge);

  cells = [pos randn(curr_n, opts.nprops)*opts.cell_variation + opts.cell_size];

  return;
end
