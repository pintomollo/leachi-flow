function find_zooids(img, systems)

  if (size(img, 3) > 1)
    img = rgb2gray(img);
  end

  perp = perpendicular_sampling(img, systems)

  if (~iscell(perp))
    perp = {perp};
  end

  for i=1:length(perp)
    gaussian_mex
    local_extrema
  end

  keyboard

  return;
end
