function [data, opts] = simulate_flow(opts)

  if (nargin == 0)
    opts = get_struct('simulation');
  end

  store_data = false;
  if (nargout > 0)
    store_data = true;
  end

  if (~isempty(opts.init_simulation) && isempty(opts.creation_params))
    opts = opts.init_simulation(opts);
  end

  cells = opts.create_cells(opts);
  opts = opts.move_cells(cells, 0, opts);

  nframes = ceil(opts.duration/opts.dt) + 1;
  frames = [0:nframes]*opts.dt;

  if (store_data)
    data = absolutepath(get_new_name('flow(\d+)\.ome\.tiff?', 'TmpData'));
    hwait = waitbar(0,'Simulating flow...','Name','B. leachi blood flow');
    img_params = [];
  else
    figure;
    colormap(gray);
  end

  thresh = 0;
  for nimg = 1:nframes
    img = draw_gaussians_mex(opts.image_size, cells) + randn(opts.image_size)*opts.image_noise;

    curr_t = frames(nimg);
    c = 0;
    while (curr_t < frames(nimg+1))
      [dmov, dt] = opts.move_cells(cells, frames(nimg), opts);
      [cells, thresh] = opts.remesh_cells(cells, dmov, thresh, opts);

      curr_t = curr_t + dt;
      c = c + 1;
    end

    if (store_data)
      [img, img_params] = all2uint16(img, img_params);

      save_data(data, img);
      waitbar(nimg/nframes,hwait);
    else
      imagesc(1024-img, [0 1024]);
      if (thresh==0)
        title(['remesh ' num2str(c)])
      else
        title(num2str(c))
      end
      drawnow;
    end
  end

  if (store_data)
    close(hwait);
  end

  return;
end
