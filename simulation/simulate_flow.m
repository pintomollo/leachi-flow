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
    img = draw_cells(opts.image_size, cells, opts) + randn(opts.image_size)*opts.image_noise;

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

function new_img = draw_cells(img_size, cells, opts)

  if (opts.ndims == 2)
    new_img = zeros(img_size);

    X = repmat([1:img_size(2)], img_size(1), 1);
    Y = repmat([1:img_size(1)].', 1, img_size(2));

    isigma = -1 ./ (2 * cells(:,3).^2);

    if (opts.nprops == 1)
      for i = 1:size(cells, 1)
        new_img = new_img + exp(isigma(i) * ((X - cells(i,1)).^2 + (Y - cells(i,2)).^2));
      end
    elseif (opts.nprops == 2)
      for i = 1:size(cells, 1)
        new_img = new_img + cells(i, 4) * exp(isigma(i) * ((X - cells(i,1)).^2 + (Y - cells(i,2)).^2));
      end
    else
      error('Not implemented yet !')
    end
  else
    error('Not implemented yet !')
  end

  return;
end
