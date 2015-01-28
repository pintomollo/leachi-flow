function [mytracking, opts] = preprocess_movie(mytracking, opts)
% PREPROCESS_MOVIE converts the OME-TIFF recordings contained in a tracking structure
% into properly filtered (as defined by the structure, see inspect_channels.m) UINT16
% files.
%
%   [MYTRACKING] = PREPROCESS_MOVIE(MYTRACKING, OPTS) rescales all the recordings used
%   in the tracking experiement using OPTS.
%
%   [MYTRACKING, OPTS] = PREPROCESS_MOVIE(...) returns in addition OPTS.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 18.05.2014

  % Initialize some computer-specific variables required for the conversion
  maxuint = intmax('uint16');

  % A nice status-bar if possible
  if (opts.verbosity > 1)
    hwait = waitbar(0,'','Name','Cell Tracking','Visible','off');
  end

  % Get the number of channels to parse
  nchannels = length(mytracking.channels);

  % Loop over all of them
  for k = 1:nchannels

    % We need the absolute path for Java to work properly
    fname = absolutepath(mytracking.channels(k).fname);

    % Hide the bar in case we are loop several times
    if (opts.verbosity > 1)
      set(hwait, 'Visible','off');
    end

    % Now we extract the corresponding metadata for potential later use
    curdir = pwd;
    cmd_path = which('bfconvert.bat');

    % We need LOCI to do so...
    if (isempty(cmd_path))
      error('Tracking:lociMissing', 'The LOCI command line tools are not present !\nPlease follow the instructions provided by install_cell_tracking');
    end
    [mypath, junk] = fileparts(cmd_path);

    % This can take a while, so inform the user
    hInfo = warndlg('Populating metadata, please wait.', 'Preprocessing movie...');

    % Move to the correct folder
    cd(mypath);

    % And call the LOCI utility to extract the metadata
    if (ispc)
      cmd_name = ['"' fname '"'];
      [res, metadata] = system(['showinf.bat -nopix -nometa -omexml-only ' cmd_name]);

    else
      cmd_name = strrep(fname,' ','\ ');
      [res, metadata] = system(['./showinf -nopix -nometa -omexml-only ' cmd_name]);
    end

    % Delete the information if need be
    if (ishandle(hInfo))
      delete(hInfo);
    end

    % Go back to the original folder
    cd(curdir);

    % Check if an error occured
    if (res ~= 0)
      error(metadata);
    end

    % Store the resulting metadata
    mytracking.channels(k).metadata = metadata;

    % Stire the original file name as we will replace it by the rescaled one
    mytracking.channels(k).file = mytracking.channels(k).fname;

    % Perfom some string formatting for the display
    indx = strfind(mytracking.channels(k).file, filesep);
    if (isempty(indx))
      indx = 1;
    else
      indx = indx(end) + 1;
    end
    if (opts.verbosity > 1)
      waitbar(0, hwait, ['Preprocessing Movie ' strrep(mytracking.channels(k).file(indx:end),'_','\_')]);
      set(hwait, 'Visible', 'on');
    end

    % Get the absolute file name
    fname = absolutepath(mytracking.channels(k).file);

    % Get the name of the new file
    tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));

    % Get the number of frames
    nframes = size_data(fname);

    % Temporary parameters about the type of data contained in the reader
    img_params = [];

    % Loop over the frames
    for i=1:nframes
      % Convert the image into UINT16
      [img, img_params] = all2uint16(load_data(fname, i), img_params);

      % Perform the required filtering
      if (mytracking.channels(k).detrend)
        img = imdetrend(img, opts.filtering.detrend_meshpoints);
      end
      if (mytracking.channels(k).cosmics)
        img = imcosmics(img, opts.filtering.cosmic_rays_window_size, opts.filtering.cosmic_rays_threshold);
      end
      if (mytracking.channels(k).hot_pixels)
        img = imhotpixels(img, opts.filtering.hot_pixels_threshold);
      end

      % Get the current range of values
      minimg = min(img(:));
      maximg = max(img(:));

      % We'll store the biggest range, to rescale it afterwards
      if(minimg < mytracking.channels(k).min)
        mytracking.channels(k).min = minimg;
      end
      if(maximg > mytracking.channels(k).max)
        mytracking.channels(k).max = maximg;
      end

      % Save the image in the temporary file
      save_data(tmp_fname, img);

      % Update the progress bar if needed
      if (opts.verbosity > 1)
        if (mytracking.channels(k).normalize)
          waitbar(i/(2*nframes),hwait);
        else
          waitbar(i/nframes,hwait);
        end
      end
    end

    % Rescale if required by the user
    if (mytracking.channels(k).normalize)
      % Get a third file to write into
      fname = tmp_fname;
      tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
      mytracking.channels(k).fname = tmp_fname;

      % Loop again over the frames
      for i=1:nframes

        % Load and rescale using the previously measured range
        img = load_data(fname, i);
        img = imnorm(img, mytracking.channels(k).min, mytracking.channels(k).max, '', 0, maxuint);

        % And save the final image
        save_data(tmp_fname, img);

        % Update the progress bar
        if (opts.verbosity > 1)
          waitbar(0.5 + i/(2*nframes),hwait);
        end
      end

      % Delete the intermidary file (i.e. the filtered one)
      delete(fname);
    else
      mytracking.channels(k).fname = tmp_fname;
    end
  end

  % Close the status bar
  if (opts.verbosity > 1)
    close(hwait);
  end

  return;
end
