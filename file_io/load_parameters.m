function opts = load_parameters(opts, fnames)
% LOAD_PARAMETERS loads parameters from a configuration file into the options structure.
%
%   OPTS = LOAD_PARAMETERS(OPTS, FNAME) loads the parameters listed in FNAME into OPTS.
%   For an example of the syntax of configuration files, see Config/default_params.txt.
%
%   OPTS = LOAD_PARAMETERS(OPTS) loads the parameters from OPTS.CONFIG_FILE or prompts
%   the user to select interactively a file.
%
%   OPTS = LOAD_PARAMETERS(FNAME) uses the standard options structure provided by
%   get_struct('options').
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  % In case there is only one argument, it might either be opts or fname
  if (nargin == 1)

    % If it contains this field, then it's opts
    if (isfield(opts, 'config_files'))
      fnames = opts.config_files(:);
      opts.config_files = {};

    % If it's a string, then it's fname
    elseif (ischar(opts) || iscell(opts))
      fnames = opts;
      opts = get_struct('options');

    % Otherwise, we just don't do anything
    else
      return;
    end
  end

  % If the name is not provided, we try to ask for one
  if (isempty(fnames))

    % Check if we can save them in the Config folder
    if (exist('Config', 'dir'))
      conf_dir = which('Config.');
    elseif (exist(['leachi-flow' filesep 'Config'], 'dir'))
      conf_dir = ['leachi-flow' filesep 'Config'];
    else
      conf_dir = pwd;
    end

    % Ask the user for a filename
    [fnames, pathname] = uigetfile({'*.txt','All text files'; '*.*','All files' }, ...
                                   'Load parameters', [conf_dir filesep]);

    % This means the user has canceled
    if (all(fnames == 0))
      return;
    end

    % Get the full name
    fnames = fullfile(pathname, fnames);
  end

  % If opts itself is some text, than we get the corresponding structure
  if (ischar(opts))
    opts = get_struct(opts);
  end

  % Format to have everything as a cell list, in case several are provided
  if (ischar(fnames))
    fnames = {fnames};
  end

  % Loop and process each file sequentially
  for i = 1:length(fnames)
    fname = fnames{i};
    fnames{i} = relativepath(fname);

    % If the file does not exists, we have a few other options
    if (~exist(fname, 'file'))

      % Maybe the extension was forgotten
      if (exist([fname '.txt'], 'file'))
        fname = [fname '.txt'];

      % Maybe it's located in the configuration folder
      elseif (exist(['Config' filesep fname], 'file'))
        fname = ['Config' filesep fname];

      % Or maybe even both previous cases
      elseif (exist(['Config' filesep fname '.txt'], 'file'))
        fname = ['Config' filesep fname '.txt'];

      % Maybe it's located in the configuration sub-folder
      elseif (exist(['leachi-flow' filesep 'Config' filesep fname], 'file'))
        fname = ['leachi-flow' filesep 'Config' filesep fname];

      % Or maybe even both previous cases
      elseif (exist(['leachi-flow' filesep 'Config' filesep fname '.txt'], 'file'))
        fname = ['leachi-flow' filesep 'Config' filesep fname '.txt'];

      % Otherwise we ran out of options
      else
        warning('Bleachi:load_parameters', ['Configuration file ''' fname ''' could not be found.'])

        continue;
      end
    end

    % We open it in text mode 
    fid = fopen(fname,'rt');

    % If there is an error, maybe we don't have the absolute path
    if (fid<0)
      fname = fullfile(pwd, fname);

      % And if it still does not work, then we skip this file
      fid = fopen(fname,'rt');
      if (fid<0)
        continue;
      end
    end

    % We can have prefixes to access subfields
    prefix = '';

    % We loop throughout the file, line by line
    line = fgetl(fid);
    while ischar(line)

      % We remove unsignificant white spaces
      line = strtrim(line);

      % We ignore empty lines
      if (length(line) > 0)

        % Prefixes start by a '#'
        if (line(1) == '#')

          % Need to avoid an error if the prefix is null
          if (length(line) > 1)

            % Extract the prefix
            prefix = line(2:end);

            % And we need a trailing dot to correctly access subfields
            if (prefix(end) ~= '.')
              prefix = [prefix '.'];
            end

          % Null prefix to reaccess fields at the root of opts
          else
            prefix = '';
          end

        % We avoid also comments which starts by '%'
        elseif (line(1) ~= '%')

          % We extract the field name (only chars) and the corresponding value
          tokens = regexp(line,'^(\S+)\s+(.+)$','tokens');

          % If we found the two elements, we can assign the value to the field
          if (~isempty(tokens) && length(tokens{1}) == 2)

            % Extract the values
            field = tokens{1}{1};
            value = tokens{1}{2};

            % Check if empty
            is_empty = strncmp(value, '[]', 2);

            % Check if an array size is provided
            tokens = regexp(value,'^(.+)#(\[[0-9 ]+\])$','tokens');
            resize = false;

            % If so, split the array and the size
            if (~isempty(tokens) && length(tokens{1}) == 2)
              value = tokens{1}{1};
              ssize = tokens{1}{2};

              resize = true;
            end

            % We use the eval function to interpret the values as in MATLAB 
            try
              if (is_empty)
                eval(['opts.' prefix field '(1:end) = [];']);
              elseif (resize)
                eval(['opts.' prefix field ' = reshape(' value ', ' ssize ');']);
              else
                eval(['opts.' prefix field ' = ' value ';']);
              end
            catch ME

              % There is an exception for empty fields
              if (is_empty)
                try
                  eval(['opts.' prefix field ' = [];']);
                catch ME
                  warning('Bleachi:load_parameters', ['An error occured when loading field ''' prefix '.' field '''\n' ME.message])
                  break;
                end
              else
                warning('Bleachi:load_parameters', ['An error occured when loading field ''' prefix '.' field '''\n' ME.message])

                break;
              end
            end
          end
        end
      end

      % Process the next line
      line = fgetl(fid);
    end
    fclose(fid);
  end

  % Finally we store the used configuration files in the provided structure, if
  % it contains the proper field
  if (isfield(opts, 'config_files'))
    if (isempty(opts.config_files))
      opts.config_files = fnames;
    elseif (ischar(opts.config_files))
      opts.config_files = [{opts.config_files}; fnames(:)];
    else
      opts.config_files = [opts.config_files(:); fnames(:)];
    end

    % And remove duplicates as well as empty calls
    files = opts.config_files(~cellfun('isempty', opts.config_files));
    [values, indx, junk] = unique(files);
    indx = sort(indx);
    opts.config_files = files(indx);
  end

  % Set the new lines back to normal
  opts = set_new_lines(opts);

  % Recompute the pixel size just in case
  opts = set_pixel_size(opts);

  return;
end

function mystruct = set_new_lines(mystruct)

  if (~isempty(mystruct))

    if (ischar(mystruct))
      mystruct = regexprep(mystruct, '\\n', '\n');
    elseif (iscell(mystruct))
      for i=1:numel(mystruct)
        mystruct{i} = set_new_lines(mystruct{i});
      end
    elseif (isstruct(mystruct))
      % In addition to the different fields, structures can be arrays.
      fields = fieldnames(mystruct);
      for i=1:numel(mystruct)

        % Loop through the fields and call myprint again, adapting the prefix
        for j=1:length(fields)
          mystruct(i).(fields{j}) = set_new_lines(mystruct(i).(fields{j}));
        end
      end
    end
  end

  return;
end
