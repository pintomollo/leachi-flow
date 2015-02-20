function clean_tmp_files()
% CLEAN_TMP_FILES removes all unused data in TmpData by recursively parsing the
% content of the .mat files present in the current directory.
%
%   [] = CLEAN_TMP_FILES() removes all unused temporary files in TmpData.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 20.12.2014

  % List the matlab files present in the current directory
  ls_dir = dir('*.mat'); 
  nfiles = length(ls_dir);

  % Prepare the list of useful temporary files
  used_tmp = [];

  % Some fancy progress status bar
  hwait = waitbar(0,'','Name','Cleaning TmpData','Visible','on');
  waitbar(0, hwait, ['Parsing .mat files to identify used data...']);

  % Loop through all the .mat files, looking for used temporary files
  for i = 1:nfiles

    % Load the current .mat file
    data = load(ls_dir(i).name);

    % Check if it contains tracking data
    if (~isfield(data, 'myrecording'))
      continue;
    end

    % Parse the whole structure and extract the ID of the used files
    tmp_tmp = recursive_tmp(data.myrecording);
    used_tmp = [used_tmp; tmp_tmp];

    % Update hte waitbar
    waitbar(i/nfiles,hwait);
  end

  % Keep only single-copied of each file
  used_tmp = unique(used_tmp);

  % Find the location of TmpData
  cd_dir = [pwd filesep];
  if (exist('TmpData', 'dir'))
    tmp_dir = [cd_dir 'TmpData' filesep];
  else
    tmp_dir = cd_dir;
  end

  % Get the content of TmpData
  ls_dir = dir([tmp_dir 'tmpmat*']);
  nfiles = length(ls_dir);

  % Update the waitbar
  waitbar(0, hwait, ['Parsing TmpData to identify existing data...']);

  % Prepare the list of files to be removed
  removed_files = {};

  % Loop over the files present in TmpData and check whether they are used
  % by any .mat file
  for d = 1:nfiles

    % Extract the ID of the corresponding temporary file
    [tmp tokens] = regexp(ls_dir(d).name, 'tmpmat(\d+)\.*', 'match', 'tokens');

    % Do we have a valid file ?
    if (length(tokens) ~= 0)

      % Get the ID and check whether it is in the list of used files
      tmp_token = str2double(char(tokens{1}));
      if (~any(used_tmp == tmp_token))

        % Otherwise, a new candidate for removal !
        removed_files(end+1) = {ls_dir(d).name};
      end
    end

    % Update the waitbar
    waitbar(d/nfiles,hwait);
  end

  % Hide it while we ask for confirmation
  set(hwait, 'Visible', 'off');

  % Do we have any files to delete ?
  if (~isempty(removed_files))

    % Ask a confirmation from the user
    answer = questdlg([{'The following files will be permanently deleted ' ...
                        'because they are not used by any project file:'}; ...
                        {''}; removed_files(:); {''}; 'Proceed ?'], ...
                        'Confirm deletion of temporary files');

    % If he agrees, delete all the listed files
    if (strcmp(answer, 'Yes'))

      % Update the waitbar
      waitbar(0, hwait, ['Deleting temporary data...']);
      set(hwait, 'Visible', 'on');

      % And loop over the files to be deleted
      nfiles = length(removed_files);
      for i = 1:nfiles

        % The actual deletion !
        delete([tmp_dir removed_files{i}]);

        % Updating the waitbar again
        waitbar(i/nfiles,hwait);
      end
    end

  % Otherwise, congratulations, it is already clean !
  else
    warndlg('TmpData is already clean !', 'Cleaning TmpData');
  end

  % Close the waitbar
  close(hwait);

  return;
end

function tmps = recursive_tmp(mystruct)

  tmps = [];
  used_tmp = [];

  fields = fieldnames(mystruct);

  indx = ismember(fields, 'fname');

  if (any(indx))
    [tmp tokens] = regexp(mystruct.fname,'tmpmat(\d+)\.*','match','tokens');
    if (length(tokens)~=0)
      used_tmp = [used_tmp; str2double(char(tokens{1}))];

      if (~exist(mystruct.fname, 'file'))
        display(['Missing ' mystruct.fname]);
      end
    end

    fields = fields(~indx);
  end

  for j=1:length(fields)
    if (isstruct(mystruct.(fields{j})))
      for k=1:length(mystruct.(fields{j}))
        if (k==1)
          sub_fields = fieldnames(mystruct.(fields{j})(k));
          has_struct = any(ismember(sub_fields, 'file'));
          if (~has_struct)
            for s = 1:length(sub_fields)
              if (isstruct(mystruct.(fields{j})(k).(sub_fields{s})))
                has_struct = true;
                break;
              end
            end
          end
        end

        if (has_struct)
          tmp_tmp = recursive_tmp(mystruct.(fields{j})(k));

          tmps = [tmps; tmp_tmp];
        end
      end
    end
  end

  tmps = [tmps; used_tmp];

  return;
end
