function [converted_file] = convert_movie(name)
% CONVERT_MOVIE converts a recording such that it can be tracked properly.
%
%   [CONVERTED] = CONVERT_MOVIE(NAME) converts the movie NAME into CONVERTED
%   which is a single OME-TIFF stack file for convenience. The LOCI command toolbox
%   is used for this process. This allows to work with a single file type later on.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 04.07.2014

  % Initialization
  converted_file = '';
  curdir = '';
  dirpath = '';
  if (nargin == 0)
    name = '';
  end

  % We got the file to be loaded as input !
  if (ischar(name) & ~isempty(name))

    % Chech whether the file is in a subdirectory
    indx = strfind(name, filesep);

    % If this is the case, separate the path and the name
    if (~isempty(indx))
      dirpath = name(1:indx(end));
      fname = name(indx(end)+1:end);

      % Check whether this is an absolute path or not
      if (dirpath(1) ~= filesep && isempty(strfind(name, ':')))
        dirpath = ['.' filesep dirpath];
      end

      % Reconstruct the full path
      fname = fullfile(dirpath, fname);
    else
      fname = name;
    end
  else

    % In case a subfolder name Movies exists, move into it for prompting
    curdir = '';
    if(exist('Movies', 'dir'))
      curdir = cd;
      cd('Movies');
    elseif(exist(['..' filesep 'Movies'], 'dir'))
      curdir = cd;
      cd(['..' filesep 'Movies']);
    end

    % Fancy output
    disp('[Select a movie file]');

    % Prompting the user for the movie file
    [fname, dirpath] = uigetfile({'*.*'}, ['Load a movie']);
    fname = fullfile(dirpath, fname);
  end

  % Return back to our original folder
  if(~isempty(curdir))
    cd(curdir);
  end

  % If no file was selected, stop here
  if (length(fname) == 0  ||  isequal(dirpath, 0))
    disp(['No movie selected']);
    return;
  end

  % Check the type of the movie chosen
  [file_path, filename, ext] = fileparts(fname);

  % If the user chose a MAT file, load it and stop here
  if (strncmp(ext, '.mat', 4))
    load(fname);
    return;
  end

  % We convert the provided type into a more handy one using LOCI
  converted_file = bftools_convert(fname);

  return;
end

% Here we'll call the LOCI Bio-formats toolbox to convert anything into our
% favorit OME-TIFF format
function [newfile] = bftools_convert(fname)

  % We need the absolute path for Java to work properly
  fname = absolutepath(fname);

  if (exist(fname, 'file') ~= 2)
    error('Tracking:BadFile', ['File ' fname ' does not exist']);
  end

  % Split the filename
  [file_path, filename, ext] = fileparts(fname);

  % Remove the extension
  name = fullfile(file_path, filename);

  % Identify the filename VS the path
  [slash] = findstr(name, filesep);
  if(length(slash)>0)
    slash = slash(end) + 1;
  else
    slash = 1;
  end

  % Creat the fancy name for display (otherwise it thinks they are LaTeX commands)
  printname = strrep(name(slash:end),'_','\_');

  % Look for the LOCI command line tool
  curdir = pwd;
  cmd_path = which('bfconvert.bat');
  if (isempty(cmd_path))
    error('Tracking:lociMissing', 'The LOCI command line tools are not present !\nPlease follow the instructions provided by install_cell_tracking');
  end
  [mypath, junk] = fileparts(cmd_path);

  % This can take a while, so inform the user
  hInfo = warndlg('Parsing metadata, please wait.', 'Converting movie...');

  % Move to the correct folder
  cd(mypath);

  % And call the LOCI utility to extract the metadata
  if (ispc)
    cmd_name = ['"' fname '"'];
    [res, metadata] = system(['showinf.bat -stitch -nopix -nometa ' cmd_name]);
  else
    cmd_name = strrep(fname,' ','\ ');
    [res, metadata] = system(['./showinf -stitch -nopix -nometa ' cmd_name]);
  end

  % Delete the information if need be
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Check if an error occured
  if (res ~= 0)
    cd(curdir);
    error(metadata);
  end

  % Extract the three important informations from the extracted metadata
  format = regexp(metadata, 'file format\s*\[([ -\w]+)\]', 'tokens');
  is_rgb = regexp(metadata, 'RGB\s*=\s*(\w+)', 'tokens');
  file_pattern = regexp(metadata, 'File pattern = ([^\n]*)\n', 'tokens');

  % In case of multiple files, regroup them into one single file
  merge_cmd = '-stitch ';
  do_merge = false;
  if (~isempty(file_pattern))
    orig_pattern = file_pattern{1};
    file_pattern = regexprep(orig_pattern, '<\d+-\d+>', '');

    % In case we did not delete the pattern, it means there was nothing to merge !
    if (length(orig_pattern{1}) == length(file_pattern{1}))
      file_pattern = '';
    else

      % Otherwise, make sure they should be merged !!
      answer = questdlg(['There are several files with the naming pattern: ''' orig_pattern{1} '''. Should we merge them together ?'], 'Merging multiple files ?');
      do_merge = strncmp(answer,'Yes',3);
      if (~do_merge)
        merge_cmd = '';
        file_pattern = '';
      end
    end
  end

  % Something went terribly wrong...
  if (isempty(format) | isempty(is_rgb))
    cd(curdir);
    error('Tracking:lociFormat', ['The metadata does not present the expected information: ''file format'' and ''RGB'' :\n\n' metadata]);
  end

  % Get the information out of the search results
  format = format{1}{1};
  is_rgb = strncmp(is_rgb{1}{1}, 'true', 4);

  if (strncmpi(format,'OME-TIFF',8) && isempty(file_pattern))

    % If it's already an OME-TIFF file, we do not need to convert it, so stop here
    newfile = fname;
    cd(curdir);
    newfile = relativepath(newfile);

    return;
  end

  % RGB will not work
  if (is_rgb)
    warning('Tracking:RGB','RGB channels will be separated, please make sure that no information is lost !')
  end

  % We create an OME-TIFF file
  if (isempty(file_pattern))
    newname = [name '.ome.tiff'];
  else
    [file_path, file_name, file_ext] = fileparts(file_pattern{1});

    % If we merge the files, use also the folder name
    if (do_merge)
      [junk, tmp_name, junk] = fileparts(file_path);
      file_name = [tmp_name '_' file_name];
    end
    newname = fullfile(file_path, [file_name '.ome.tiff']);
  end

  % If the file already exists, we ask what to do
  if(exist(newname,'file'))

    % We initially do not know what to do
    answer = 0;

    % We do not accept "empty" answers
    while (answer == 0)
      answer = menu(['The OME-TIFF version of ' printname ' already exists, overwrite it ?'],'Yes','No');
    end

    % Act accorindly to the policy
    switch answer

      % Delete the current files (did not dare overwriting it directly)
      case 1
        delete(newname);

      % Otherwise we can stop here
      case 2
        % Store the new name
        newfile = newname;
        cd(curdir);
        newfile = relativepath(newfile);

        return;
    end
  end

  % This also takes quite some time, so warn the user
  hInfo = warndlg('Converting to OME-TIFF, please wait.', 'Converting movie...');

  % Call directly the command line tool to do the job
  if (ispc)
    cmd_newname = ['"' newname '"'];
    [res, infos] = system(['bfconvert.bat ' merge_cmd '-separate ' cmd_name ' ' cmd_newname]);
  else
    cmd_newname = strrep(newname,' ','\ ');
    [res, infos] = system(['./bfconvert ' merge_cmd '-separate ' cmd_name ' ' cmd_newname]);
  end

  % Delete the window if need be
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Check if an error occured
  if (res ~= 0)
    cd(curdir);
    error(infos);
  end

  % Store the new name in relative path and come back to the original folder
  newfile = newname;
  cd(curdir);
  newfile = relativepath(newfile);

  return;
end
