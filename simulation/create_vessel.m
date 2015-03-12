function opts = create_vessel(opts)

  b_leachi = get_struct('botrylloides_leachi');
  vessel_props = gmdistribution(b_leachi.vessel_width.mu, b_leachi.vessel_width.sigma, b_leachi.vessel_width.proportions);
  vessel_widths = random(vessel_props, opts.init_params(1));

  vessels = NaN(opts.init_params(1), 3);
  for i=1:opts.init_params(1)
    vessels(i, :) = [rand(1)*(2*pi) vessel_widths(i) rand(1)*opts.image_size(1)];
  end

  opts.creation_params = vessels(1, 1:2);
  opts.movement_params = [cos(vessels(1,1)) sin(vessels(1, 1))];

  return;
end
