function data = simulate_flow(opts)

  if (nargin == 0)
    opts = get_struct('simulation');
  end

  store_data = false;
  if (nargout > 0)
    store_data = true;
  end

  cells = opts.create_cells(opts);

  frames = [0:opts.dt:opts.duration];
  nframes = length(frames);

  if (store_data)
    data = zeros([opts.image_size nframes]);
  else
    figure;
  end

  for nimg = 1:nframes
    img = draw_cells(opts.image_size, cells) + randn(opts.image_size)*opts.image_noise;

    cells = opts.move_cells(cells, frames(nimg), opts);

    if (store_data)
      data(:,:,nimg) = img;
    else
      imagesc(img);
      drawnow;
    end
  end

  return;
end

function new_img = draw_cells(img_size, cells)

  new_img = zeros(img_size);

  X = repmat([1:img_size(2)], img_size(1), 1);
  Y = repmat([1:img_size(1)].', 1, img_size(2));

  isigma = -1 ./ (2 * cells(:,3).^2);

  for i = 1:size(cells, 1)
    new_img = new_img + exp(isigma(i) * ((X - cells(i,1)).^2 + (Y - cells(i,2)).^2));
  end

  return;
end
