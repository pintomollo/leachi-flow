function parameters = For3D(varargin)

  % Inputs processing
  if (length(varargin) > 0 && isstruct(varargin{1}))
    parameters = update_structure(varargin{1}, 'For3D');
    varargin(1) = [];
  else
    parameters = get_struct('For3D');
  end

  % Now we check that the parameters were provided in pairs
  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error 'Properties pairs must come in PAIRS.';
  end

  % Loop over the pairs of parameters
  for i = 1:npairs
    % If the parameter exists in opts we simply assign it the
    % provided value
    if (isfield(parameters, varargin{2*i - 1}))
      parameters.(varargin{2*i - 1}) = varargin{2*i};

    % Or ignore it !
    else
      warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
    end
  end

  % First make sure the parameters are correct
  parameters = get_parameters(parameters);

  % Resample the images first !
  parameters = resample_tif(parameters);

  % Then adjust the image into a proepr stack
  parameters.filename = adjust_tif(parameters.filename);
  parameters.filename = resize_tif(parameters.filename);

  % The splitting of colors stuff
  if (parameters.colorize)
    parameters.filename = colorize_stack(parameters.filename);
  end

  % Then align the stack
  parameters.filename = register_stack(parameters.filename, parameters.min_fraction);

  % Filter the stack
  parameters.filename = smooth_slices(parameters.filename);
  parameters.filename = equalize_stack(parameters.filename, parameters.alpha);
  parameters = filter_stack(parameters);

  % Reconstruct the volume
  %rendering_3D(parameters);

  if (nargout == 0)
    clear parameters;
  end

  return;
end

% Go through the parameters
function params = get_parameters(params)

  % Get the files provided by the user
  files = get_filenames(params.filename);

  % If none, ask for some
  N = length(files);
  if (N == 0)
    [fname, pathname] = uigetfile('*.*', 'Select one of the section that you want to reconstruct.');

    if isequal(fname,0)
      error('User canceled the selection.');
    else
      [filepath, fname, fileext] = fileparts(fname);
      new_name = fullfile(pathname, fullfile(filepath, ['*' fileext]));

      files = get_filenames(new_name);
      N = length(files);

      % Keep them only if there are some
      if (N > 0)
        params.filename = new_name;
      else
        error('No valid image selected.');
      end
    end
  end

  % Try reading the metadata
  info = imfinfo(files{1});

  % For now on, I look only for one field
  if (isfield(info, 'ImageDescription'))
    xval = regexp(info.ImageDescription, 'XCalibrationMicrons=([\d\.]+)', 'tokens');

    if (~isempty(xval))
      yval = regexp(info.ImageDescription, 'YCalibrationMicrons=([\d\.]+)', 'tokens');

      if (~isempty(yval))
        resol = 0.5*(str2double(xval{1}) + str2double(yval{1}));

        if (isfinite(resol))
          params.pixel_size = resol;
        end
      end
    end
  end

  % Ask a confirmation from the user
  params = edit_options(params);
  drawnow;

  return;
end
