function params = resample_tif(params)

  if nargin < 1
    params = get_struct('For3D');
  end

  files = params.filename;

  new_names = {};
  dir_out = '_resampled'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  center = ceil(N/2);

  Nsampling = (params.slice_width/params.pixel_size/1.5); % assuming ideally dz = 1.5dx for resoluting voxel after sampling.
  if Nsampling <= 1
    nframes = ceil(1/Nsampling);
    indx = [fliplr([center:-nframes:1]) center+nframes:nframes:N];

    files = files(indx);
    params.slice_width = params.slice_width*nframes;

    N = length(files);
    center = ceil(N/2);
  end % Thymus (screen capture at 2x): floor(22/6.5/1.5) = 2, Thymus (3x): floor(22/4.33/1.5) = 3, LN (8x): floor(22/2.6/1.5) = 5

  Nsampling = ceil(params.Nsampling);

  if (Nsampling > 1)
    indx = [fliplr([center:-Nsampling:1]) center+Nsampling:Nsampling:N];
  else
    params.Nsampling = 1;
    indx = [1:N];
  end

  resampling = params.Nsampling*(params.slice_width/params.pixel_size/1.5); % assuming ideally dz = 1.5dx for resoluting voxel after sampling.
  invsample = 1/resampling;

  fprintf(' Resampling images :     ');

  files = files(indx);
  N = length(files);
  new_names = files;

  for i = 1:N % loop over images to resize images
    fprintf('\b\b\b%3d',i);

    filename = files{i};
    im = imread(filename);

    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    im = imresize(im, invsample);

    imwrite(im, new_name, 'Compression', 'none');

    new_names{i} = new_name;
  end
  fprintf(1, '\b\b\b\bdone\n');

  params.filename = new_names;
  params.pixel_size = params.pixel_size * resampling;

  return;
end
