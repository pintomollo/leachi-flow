function vessel = create_background(vessel, simul, opts)

  if (length(simul.init_params) > 2)
    b_leachi = get_struct('botrylloides_leachi');

    background = get_struct('background', 1);
    background.vessel = vessel.border;
    background.vessel_properties = [min(b_leachi.blood_cell.mu(:,2))/2 1.5/opts.pixel_size];

    if (simul.init_params(3) > 0)
      simul.cell_density = simul.cell_density * simul.init_params(3);
      cells = simul.create_cells(simul, opts);

      cell_props = gmdistribution(b_leachi.blood_clot.mu, b_leachi.blood_clot.sigma, b_leachi.blood_clot.proportions);
      cells = [cells(:,1:2) random(cell_props, size(cells, 1))];
      cells = cells(~any(cells(:,3:4)<=0, 2), :);
      cells(:,3) = cells(:,3) / opts.pixel_size;

      background.clot = cells;
    end

    if (length(simul.init_params) > 3 && simul.init_params(4) > 0)
      cell_props = gmdistribution(b_leachi.ampulla.mu, b_leachi.ampulla.sigma, b_leachi.ampulla.proportions);

      curr_size = (simul.image_size-1) * (1 + 2*simul.outside_ridge);
      curr_n = ceil(simul.init_params(4));

      pos = bsxfun(@minus, bsxfun(@times, rand(curr_n, simul.ndims), curr_size) + 1, ...
                           (simul.image_size-1)*simul.outside_ridge);

      cells = [pos random(cell_props, curr_n)];
      cells = cells(~any(cells(:,3:4)<=0, 2), :);

      background.ampulla = cells;
    end

    vessel.background = background;
  end

  return;
end
