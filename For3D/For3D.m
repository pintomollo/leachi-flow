function parameters = For3D(varargin)

  continued = false;
  fprintf('For3D pipeline :\n')

  % Inputs processing
  if (length(varargin) > 0 && isstruct(varargin{1}))
    parameters = update_structure(varargin{1}, 'For3D');
    varargin(1) = [];
    continued = true;
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
  ndone = length(parameters.file_log);

  % Get the working directory
  [workdir, junk, junk] = fileparts(parameters.filename{1});
  if (continued)
    [workdir, junk, junk] = fileparts(workdir);
  end
  savename = fullfile(workdir, 'myFor3D.mat');

  % Saving the progress
  if (~continued || ndone<1)
    parameters.file_log{end+1} = parameters.filename;
    save(savename, 'parameters');
  end

  % Then adjust the image into a proper stack
  if (~continued || ndone<2)
    parameters.filename = resize_tif(parameters.filename);
    parameters.file_log{end+1} = parameters.filename;
  end
  if (~continued || ndone<3)
    parameters.filename = clean_borders(parameters.filename);
    parameters.file_log{end+1} = parameters.filename;

    % Saving the progress
    save(savename, 'parameters');
  end

  % Then align the stack
  if (~continued || ndone<4)
    parameters.filename = register_stack(parameters.filename, parameters.min_fraction);

    % Saving the progress
    parameters.file_log{end+1} = parameters.filename;
    save(savename, 'parameters');
  end

  % Smooth the volume and the intensity of the organ
  if (~continued || ndone<5)
    parameters.filename = smooth_slices(parameters.filename, parameters.smoothing_span);
    parameters.file_log{end+1} = parameters.filename;
  end
  if (~continued || ndone<6)
    parameters.filename = equalize_stack(parameters.filename, parameters.alpha);
    parameters.file_log{end+1} = parameters.filename;

    % Saving the progress
    save(savename, 'parameters');
  end

  % Then resample the images to the proper resolution
  if (~continued || ndone<7)
    parameters = resample_tif(parameters);

    % Saving the progress
    parameters.file_log{end+1} = parameters.filename;
    save(savename, 'parameters');
  end

  % Adjust them to fit an expected sparse matrix
  if (~continued || ndone<8)
    parameters.filename = adjust_tif(parameters.filename);
    parameters.file_log{end+1} = parameters.filename;

    % Saving the progress
    save(savename, 'parameters');
  end

  % The splitting of colors stuff
  if (parameters.colorize && (~continued || (ndone<8+parameters.colorize)))
    parameters.filename = colorize_stack(parameters.filename);
    parameters.file_log{end+1} = parameters.filename;

    % Saving the progress
    save(savename, 'parameters');
  end

  % Filter the stack
  if (~continued || (ndone<9+parameters.colorize))
    parameters = filter_stack(parameters);
    parameters.file_log{end+1} = parameters.filename;

    % Saving the progress
    save(savename, 'parameters');
  end

  % Reconstruct the volume
  rendering_3D(parameters);

  if (nargout == 0)
    clear parameters;
  end

  return;
end

% Go through the parameters
function params = get_parameters(params)

  % Get the files provided by the user
  files = get_filenames(params.filename);
  params.filename = files;

  % If none, ask for some
  N = length(files);
  if (N == 0)
    fprintf(' Select an image from the sections to reconstruct in 3D : ');

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
        params.filename = files;
      else
        error('No valid image selected.');
      end
    end

    fprintf('thanks !\n');
  end

  % Try reading the metadata
  info = imfinfo(files{1});

  % For now on, I look only for one field
  if (isfield(info, 'ImageDescription'))
    xval = regexp(info(1).ImageDescription, 'XCalibrationMicrons=([\d\.]+)', 'tokens');

    if (~isempty(xval))
      yval = regexp(info(1).ImageDescription, 'YCalibrationMicrons=([\d\.]+)', 'tokens');

      if (~isempty(yval))
        resol = 0.5*(str2double(xval{1}) + str2double(yval{1}));

        if (isfinite(resol))
          params.pixel_size = resol;
        end
      end
    end
  end

  % Ask a confirmation from the user
  fprintf(' Confirm the parameter values of the experiment : ');
  [params, updated] = edit_options(params);
  if (updated)
    fprintf('thanks !\n');
    drawnow;
  else
    error('Pipeline cancelled by the user.');
  end

  return;
end
