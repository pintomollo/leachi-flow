function params = resample_tif(params)

  if nargin < 1
    params = get_struct('For3D');
  end

  files = params.filename;

  Nsampling = floor(params.slice_width/params.pixel_size/1.5); % assuming ideally dz = 1.5dx for resoluting voxel after sampling.
  if Nsampling <= 1
    return;
  end % Thymus (screen capture at 2x): floor(22/6.5/1.5) = 2, Thymus (3x): floor(22/4.33/1.5) = 3, LN (8x): floor(22/2.6/1.5) = 5

  insample = 1/Nsampling;

  new_names = {};
  dir_out = '_resampled'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Resampling images...');

  new_names = files;

  for i = 1:N % loop over images to resize images
    filename = files{i};
    im = imread(filename);

    im = imresize(im, insample);

    imwrite(im, new_name, 'Compression', 'none');

    new_names{i} = new_name;
  end

  params.filename = new_names;
  params.pixel_size = params.pixel_size * Nsampling;

  return;
end
