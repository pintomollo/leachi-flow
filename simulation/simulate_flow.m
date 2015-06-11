function [data, opts] = simulate_flow(opts)

  % Get the parameter structure if not provided
  if (nargin == 0)
    opts = get_struct('simulation');
  end

  % Decide whether we display or store
  store_data = false;
  if (nargout > 0)
    store_data = true;
  end

  % Initialize the simulation
  if (~isempty(opts.init_simulation) && isempty(opts.creation_params))
    opts = opts.init_simulation(opts);
  end

  % Need the row/colum size of the image as well (most of the code uses X/Y data)
  img_size = opts.image_size([2 1]);

  % Create the cells and initialize the moving step
  cells = opts.create_cells(opts);
  opts = opts.move_cells(cells, 0, opts);

  % Get the frames
  nframes = ceil(opts.duration/opts.dt) + 1;
  frames = [0:nframes]*opts.dt;

  % Get ready to output our simulation
  if (store_data)
    data = absolutepath(get_new_name('flow(\d+)\.ome\.tiff?', 'TmpData'));
    infos = get_struct('image_infos');
    infos.scaling = double(intmax('uint16'));
    hwait = waitbar(0,'Simulating flow...','Name','B. leachi blood flow');
  else
    figure;
    colormap(gray);
  end

  % Now iterate over all the frames
  thresh = 0;
  curr_t = frames(1);
  for nimg = 1:nframes

    % Draw the current flow
    img = draw_gaussians_mex(img_size, cells) + randn(img_size)*opts.image_noise;
    img = 1 - img;

    % Perform simulations steps until we reach the next frame we need to output
    c = 0;
    did_remesh = false;
    while (curr_t < frames(nimg+1))

      % Move one step forward and remesh if necessary
      [dmov, dt] = opts.move_cells(cells, curr_t, opts);
      [cells, thresh] = opts.remesh_cells(cells, dmov, thresh, opts);

      % Update the time and counters
      curr_t = curr_t + dt;
      c = c + 1;
      did_remesh = (did_remesh || thresh==0);
    end

    % Either store or display
    if (store_data)
      img = all2uint16(img, infos);

      save_data(data, img);
      waitbar(nimg/nframes,hwait);
    else
      imagesc(img, [0 1]);
      if (did_remesh)
        title(['remesh ' num2str(c)])
      else
        title(num2str(c))
      end
      drawnow;
    end
  end

  % Might need to close the progress bar
  if (store_data)
    close(hwait);
  end

  return;
end
