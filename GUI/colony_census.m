function [] = colony_census(fname)
% INSPECT_RECORDING displays a pop-up window for the user to manually identify the
% type of data contained in the different channels of a movie recording.
%
%   INSPECT_RECORDING(CHANNELS) displays the window using the data
%   contained in CHANNELS, updates it accordingly to the user's choice and saves
%   the adequate structure for later analysis MYRECORDING. CHANNELS can either
%   be a string, a cell list of strings or a 'channel' structure (see get_struct.m).
%   MYRECORDING is a structure as defined by get_struct('myrecording').
%
%   INSPECT_RECORDING() prompts the user to select a recording and converts
%   it before opening the GUI.
%
% Blanchoud group, UNIFR
% Simon Blanchoud
% 19.12.2018

  % Argument checking, need to know if we ask for a recording or not.
  if (nargin == 0 || isempty(fname) || ...
     (isstruct(fname) && isfield(fname, 'channels') && isempty(fname.channels)))
    fname = load_images();

    % We did not get anything to handle...
    if isempty(fname)
      return;
    end
  end

  % Create the channels structure if it was not provided.
  if (isstruct(fname))
    if (isfield(fname, 'experiment'))
      myrecording = fname;
      channels = myrecording.channels;
    else
      myrecording = get_struct('myrecording');
      channels = fname;
    end
  else
    % Put everything in a cell list
    if (ischar(fname))
      fname = {fname};
    end

    % Parse the list and copy its content into the channels structure
    nchannels = length(fname);
    channels = get_struct('census', [nchannels 1]);
    for i=1:nchannels
      channels(i).fname = fname{i};
    end

    myrecording = get_struct('myrecording');
  end

  % Create the GUI
  [hFig, handles] = create_figure(channels);

  % Store all the useful data in a single structure to be stored in the figure
  data = struct('channels', channels, ...
                'recording', myrecording, ...
                'img', [], ...
                'orig_img', [], ...
                'img_next', [], ...
                'handles', handles);

  % Update its content
  data = update_display(data);

  % Store the data structure in the main figure to be retrieved by the callbacks
  guidata(hFig, data);

  % Display the figure
  set(hFig,'Visible', 'on');

  return;
end

function [hFig, handles] = create_figure(channels)
% This function actually creates the GUI, placing all the elements
% and linking the callbacks.

  % The number of channels provided
  nchannels = length(channels);

  % Initialize the possible types and compressions
  typestring = {'luminescence';'brightfield'; 'dic'; 'fluorescence'};
  typecompress = {'none', 'lzw', 'deflate', 'jpeg'};

  % Initialize the structure used for the interface
  liststring = '';
  for i = 1:nchannels

    % Build the displayed list
    liststring = [liststring 'Image ' num2str(i)];
    if (i < nchannels)
      liststring = [liststring '|'];
    end
  end

  % Create a name for the experiment based on the filename
  exp_name = channels(1).fname;
  [junk, exp_name, junk] = fileparts(exp_name);
  [junk, exp_name, junk] = fileparts(exp_name);
  exp_name = regexprep(exp_name, ' ', '');

  % Create my own grayscale map for the image display
  mygray = [0:255]' / 255;
  mygray = [mygray mygray mygray];

  % We build a list of all buttons to easily block and release them
  enabled = [];

  % Retrieve the resolution of the screen and infer the position of the GUI
  resol = get(0, 'ScreenSize');
  img_pos = ceil([0.05*resol(3:4) 0.9*resol(3:4)]);

  % The main figure, cannot be rescaled, closed nor deleted
  hFig = figure('PaperUnits', 'centimeters',  ...
                'CloseRequestFcn', @cancel_CloseRequestFcn, ...
                'Colormap', mygray, ...
                'Color',  [0.7 0.7 0.7], ...
                'MenuBar', 'none',  ...
                'Name', 'Colony census',  ...
                'NumberTitle', 'off',  ...
                'Units', 'pixels', ...
                'Position', img_pos, ...
                'DeleteFcn', @gui_Callback, ...
                'HandleVisibility', 'callback',  ...
                'Tag', 'channel_fig', ...
                'UserData', [], ...
                'Visible', 'off');

  %%%%%% Now the buttons around the main panel

  % The list of channels
  hChannel = uicontrol('Parent', hFig, ...
                       'Units', 'normalized',  ...
                       'Callback', @gui_Callback, ...
                       'Position', [0.01 0.11 0.1 0.79], ...
                       'String', liststring, ...
                       'Style', 'listbox',  ...
                       'Value', 1, ...
                       'Tag', 'channels');
  enabled = [enabled hChannel];

  % The OK button
  hOK = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @channel_CloseRequestFcn, ...
                  'Position', [0.70 0.02 0.18 0.05], ...
                  'String', 'OK',  ...
                  'Tag', 'pushbutton11');
  enabled = [enabled hOK];

  % The Cancel button
  hCancel = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @cancel_CloseRequestFcn, ...
                  'Position', [0.90 0.02 0.08 0.05], ...
                  'String', 'Cancel',  ...
                  'Tag', 'pushbutton12');
  enabled = [enabled hCancel];

  % The Add and Remove buttons
  hAdd = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @add_channel_Callback, ...
                  'Position', [0.01 0.055 0.1 0.04], ...
                  'String', 'Add image',  ...
                  'Tag', 'pushbutton13');
  enabled = [enabled hAdd];

  hRemove = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @remove_channel_Callback, ...
                  'Position', [0.01 0.01 0.1 0.04], ...
                  'String', 'Remove image',  ...
                  'Tag', 'pushbutton14');
  enabled = [enabled hRemove];

% Dragzoom help message
imghelp = regexp(help('dragzoom2D'), ...
           '([ ]+Interactions:.*\S)\s+Adapted from','tokens');
imghelp = ['DRAGZOOM2D interactions (help dragzoom2D):\n\n', imghelp{1}{1}];

  % The experiment name and its labels
  hText = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.2 0.93 0.09 0.025], ...
                    'String', 'Experiment name:',  ...
                    'TooltipString', sprintf(imghelp), ...
                    'BackgroundColor', get(hFig, 'Color'), ...
                    'FontSize', 12, ...
                    'Style', 'text',  ...
                    'Tag', 'text1');

  hName = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.3 0.93 0.5 0.05], ...
                    'String', exp_name,  ...
                    'FontSize', 12, ...
                    'Style', 'edit',  ...
                    'Tag', 'experiment');
  enabled = [enabled hName];

  %%%%%%% Now the main panel

  % The panel itsel
  hPanel = uipanel('Parent', hFig, ...
                   'Title', 'Image 1',  ...
                   'Tag', 'uipanel',  ...
                   'Clipping', 'on',  ...
                   'Position', [0.12 0.11 0.87 0.8]);

  % The two axes
  hAxes = axes('Parent', hPanel, ...
               'Position', [0.01 0.1 0.43 0.9], ...
               'DataAspectRatio', [1 1 1], ...
               'Visible', 'off',  ...
               'Tag', 'axes');

  hAxesNext = axes('Parent', hPanel, ...
               'Position', [0.45 0.1 0.43 0.9], ...
               'DataAspectRatio', [1 1 1], ...
               'Visible', 'off',  ...
               'Tag', 'axes');

  % The two radio button groups that handle which image to display
  % For the choices to be mutually exclusive, one has to put them inside
  % such uibuttongroup.

  % The Add and Remove buttons
  hControl = uicontrol('Parent', hPanel, ...
                  'Units', 'normalized',  ...
                  'Callback', @add_ROI_Callback, ...
                  'Position', [0.125 0.01 0.2 0.075], ...
                  'String', 'Add ROI',  ...
                  'Tag', 'addroi');
  enabled = [enabled hControl];

  % The Add and Remove buttons
  hControl = uicontrol('Parent', hPanel, ...
                  'Units', 'normalized',  ...
                  'Callback', @find_zooids_Callback, ...
                  'Position', [0.565 0.01 0.2 0.075], ...
                  'String', 'Find zooids',  ...
                  'Tag', 'addroi');
  enabled = [enabled hControl];

  % The type, color and compression of the channel, along with its labels
  hText = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.9 0.925 0.05 0.05], ...
                    'String', 'Channel:',  ...
                    'FontSize', 12, ...
                    'FontWeight', 'bold', ...
                    'Style', 'text',  ...
                    'Tag', 'text16');

  hText = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.885 0.875 0.075 0.05], ...
                    'String', 'Pixel size',  ...
                    'TooltipString', 'Pixel resolution in um', ...
                    'Style', 'text',  ...
                    'Tag', 'text17');

  hResol = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.89 0.845 0.07 0.04], ...
                    'String', '2.5', ...
                    'Style', 'edit',  ...
                    'Tag', 'data');
  enabled = [enabled hResol];

  hText = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.9 0.77 0.05 0.05], ...
                    'String', 'Amplitude',  ...
                    'TooltipString', 'Intensity threshold required between two zooids (-1 = automatic threshold)', ...
                    'Style', 'text',  ...
                    'Tag', 'text17');

  hAmpl = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.89 0.74 0.07 0.04], ...
                    'Style', 'edit',  ...
                    'String', '-1', ...
                    'Tag', 'data');
  enabled = [enabled hAmpl];

  % The various filters, along with their labels
  hText = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Position', [0.9 0.425 0.05 0.05], ...
                    'String', 'Filters:',  ...
                    'FontSize', 12, ...
                    'FontWeight', 'bold', ...
                    'Style', 'text',  ...
                    'Tag', 'text16');

  hNorm = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.35 0.1 0.05], ...
                         'String', 'Normalize',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'normalize');
  enabled = [enabled hNorm];

  % The Snapshot button
  hSnapshot = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @options_Callback, ...
                  'Position', [0.01 0.93 0.05 0.05], ...
                  'String', 'Snapshot',  ...
                  'Tag', 'snapshot');
  enabled = [enabled hSnapshot];

  % We store all the useful handles into a structure to easily retrieve them,
  % along with some indexes
  handles = struct('uipanel', hPanel, ...
                   'list', hChannel, ...
                   'normalize', hNorm, ...
                   'axes', [hAxes hAxesNext], ...
                   'experiment', hName, ...
                   'all_buttons', enabled, ...
                   'resolution', hResol, ...
                   'amplitude', hAmpl, ...
                   'hFig', hFig, ...
                   'img', -1, ...
                   'scatter', -1, ...
                   'roi', {{}}, ...
                   'display', [1 1], ...
                   'prev_channel', -1, ...
                   'current', 1);

  % Link both axes to keep the same information on both sides
  linkaxes(handles.axes);

  return;
end

function [data] = handle_rois(data, prev_indx, next_indx)

  handles = data.handles;

  % Extract the ROIs
  nrois = length(handles.roi);
  if (prev_indx > 0)
    systems = NaN(0, 2);
    for i=1:nrois
      systems = [systems; handles.roi{i}.getPosition(); NaN(1,2)];
    end
    data.channels(prev_indx).system = systems;
  end

  if (next_indx > 0)
    systems = data.channels(next_indx).system;
    sys = find(all(isnan(systems), 2));
    nsys = length(sys);

    prev = 1;
    for i=1:nsys
      curr_sys = systems(prev:sys(i)-1,:);
      if (i > nrois)
        handles.roi{i} = impoly(handles.axes(1), curr_sys, 'Closed', false);
      else
        handles.roi{i}.setPosition(curr_sys);
      end
      prev = sys(i)+1;
    end

    for i=nrois:-1:nsys+1
      handles.roi{i}.delete();
      handles.roi(i) = [];
    end
  end

  data.handles = handles;

  return;
end

function data = update_display(data, recompute)
% The main figure of the GUI, the one responsible for the proper display
% of its content.

  % By default we recompute everything
  if (nargin < 2)
    recompute = true;
  end

  % Retrive the stored data
  handles = data.handles;
  hFig = handles.hFig;

  % Get the indexes of the current frame and channel
  indx = handles.current;

  % Stop if no data at all
  if (indx < 1)
    return;
  end

  % Here we recompute all the filtering of the frame
  if (recompute || indx ~= handles.prev_channel)
    % Because it takes long, display it and block the GUI
    set(hFig, 'Name', 'Colony census (Filtering...)');
    set(handles.all_buttons, 'Enable', 'off');
    drawnow;
    refresh(hFig);
  end

  % If we have changed channel, we need to update the display of the buttons
  if (indx ~= handles.prev_channel)
    if (handles.prev_channel > 0)
      data.channels(handles.prev_channel).pixel_size = str2double(get(handles.resolution, 'String'));
      data.channels(handles.prev_channel).amplitude = str2double(get(handles.amplitude, 'String'));
    end

    % Set the name of the current panel
    set(handles.uipanel,'Title', ['Image ' num2str(indx)]);

    set(handles.normalize,'Value', data.channels(indx).normalize);

    set(handles.resolution, 'String', num2str(data.channels(indx).pixel_size));
    set(handles.amplitude, 'String', num2str(data.channels(indx).amplitude));

    if (numel(handles.img) > 1 && all(ishandle(handles.img)))
      data = handle_rois(data, handles.prev_channel, indx);
    end

    % And setup the indexes correctly
    handles.prev_channel = indx;
  end

  % Here we recompute all the filtering of the frame
  if (recompute)

    % Try to avoid reloading frames as much as possible
    data.orig_img = imread(data.channels(indx).fname);

    % Copy to the working variable
    data.img = data.orig_img;

    % Normalize the image ?
    if (data.channels(indx).normalize)
      data.img = imnorm(data.img);
    end
  end

  % Determine which image to display in the left panel
  img1 = data.img;
  img2 = data.img;

  % Get the location of the zooids
  zooids = data.channels(indx).zooids;
  if (isempty(zooids))
    zooids = NaN(1,2);
  end

  curr_color = 'k';

  % If we have already created the axes and the images, we can simply change their
  % content (i.e. CData)
  if (numel(handles.img) > 1 && all(ishandle(handles.img)))
    set(handles.img(1),'CData', img1);
    set(handles.img(2),'CData', img2);
    set(handles.scatter, 'XData', zooids(:,1), 'YData', zooids(:,2), 'Color', curr_color);
  else

    % Otherwise, we create the two images in their respective axes
    handles.img = image(img1,'Parent', handles.axes(1),...
                      'CDataMapping', 'scaled',...
                      'Tag', 'image');
    handles.img(2) = image(img2,'Parent', handles.axes(2), ...
                      'CDataMapping', 'scaled',...
                      'Tag', 'image');

    % Hide the axes and prevent a distortion of the image due to stretching
    set(handles.axes,'Visible', 'off',  ...
               'DataAspectRatio',  [1 1 1]);

    handles.scatter = line('XData', zooids(:,1), 'YData', zooids(:,2), 'Parent', handles.axes(2), 'Color', curr_color, 'Marker', 'o', 'LineStyle', 'none');

    % Drag and Zoom library from Evgeny Pr aka iroln
    dragzoom2D(handles.axes(1))
  end

  if (recompute || indx ~= handles.prev_channel)
    % Release the image
    set(hFig, 'Name', 'Colony census');
    set(handles.all_buttons, 'Enable', 'on');
  end

  % Update the data
  data.handles = handles;

  return
end

function remove_channel_Callback(hObject, eventdata)
% This function removes ones channel from the list

  data = guidata(hObject);
  handles = data.handles;

  % Get the current index for later
  indx = handles.current;

  % And ask for confirmation
  answer = questdlg(['Are you sure you want to remove channel ' num2str(indx) ' ?' ...
              '(No data will be deleted from the disk)'], 'Removing a channel ?');
  ok = strcmp(answer, 'Yes');

  % If it's ok, let's go
  if (ok)

    % Remove the current index and select the first one
    data.channels(indx) = [];
    handles.current = 1;

    % If it was the only one, we need to handle this
    if (isempty(data.channels))

      % Set up the indexes as empty
      handles.current = 0;
      handles.prev_channel = -1;

      % As this is long, block the GUI
      set(handles.all_buttons, 'Enable', 'off');
      set(handles.img, 'CData', []);
      set(handles.list, 'String', '');

      % And call the movie conversion function to load a new one
      new_channel = load_images();

      % Did we get something back ?
      if (~isempty(new_channel))

        % We provide basic default values for all fields
        data.channels = get_struct('census');
        data.channels.fname = new_channel;

        % We update the list of available channels
        set(handles.list, 'String', 'Image 1', 'Value', 1);
      end

    % Otherwise, delete the last traces of the deleted channel
    else
      handles.prev_channel = -1;
      tmp_list = get(handles.list, 'String');
      tmp_list(indx,:) = [];
      set(handles.list, 'String', tmp_list, 'Value', 1);
    end

  % Update the data
  data.handles = handles;

    % Release the GUI and update
    data = update_display(data, true);

  guidata(handles.hFig, data);

    set(handles.all_buttons, 'Enable', 'on');
  end

  return;
end

function add_ROI_Callback(hObject, eventdata)

  data = guidata(hObject);
  handles = data.handles;

  % As this is long, block the GUI
  set(handles.all_buttons, 'Enable', 'off');
  set(handles.hFig, 'Name', 'Colony census (Draw your ROI...)');

  h = impoly(handles.axes(1), 'Closed', false);
  if (~isempty(h))
    handles.roi{end+1} = h;
  end

  % Update the data
  data.handles = handles;
  guidata(handles.hFig, data);

  % Release the GUI
  set(handles.hFig, 'Name', 'Colony census');
  set(handles.all_buttons, 'Enable', 'on');
end

function find_zooids_Callback(hObject, eventdata)

  data = guidata(hObject);
  handles = data.handles;

  % As this is long, block the GUI
  set(handles.all_buttons, 'Enable', 'off');

  params = [8/str2double(get(handles.resolution, 'String')) str2double(get(handles.amplitude, 'String'))];

  data = handle_rois(data, handles.current, -1);

  %keyboard

  zooids = find_zooids(data.img, data.channels(handles.current).system, params);
  data.channels(handles.current).zooids = zooids;


  % Release the GUI
  data = update_display(data, false);
  guidata(handles.hFig, data);
  set(handles.all_buttons, 'Enable', 'on');
end

function add_channel_Callback(hObject, eventdata)
% This function adds a new channel to the current recording

  data = guidata(hObject);
  handles = data.handles;

  % Remember how many there were
  nchannels = length(data.channels);

  % As this is long, block the GUI
  set(handles.all_buttons, 'Enable', 'off');

  % And call the movie conversion function to get a new channel
  new_channel = load_images();

  % Did we get anything ?
  if (~isempty(new_channel))

    liststring = get(handles.list, 'String');
    for i=1:length(new_channel)
      % We provide basic default values for all fields
      data.channels(end+1) = get_struct('census');
      data.channels(end).fname = new_channel{i};

      % We update the list of available channels
      if (nchannels == 0)
        liststring = ['Image ' num2str(length(data.channels))];
      else
        if (size(liststring, 1) > 1)
          liststring(:,end+1) = '|';
          liststring = liststring.';
          liststring = liststring(:).';
          liststring = liststring(1:end-1);
        end
        liststring = [liststring '|Image ' num2str(length(data.channels))];
      end

      nchannels = nchannels + 1;
    end
    set(handles.list, 'String', liststring);
  end

  % If we just add a new first channel, we need to set it up properly and update
  if (nchannels == 0 && length(channels) > 0)
    handles.current = 1;
    set(handles.list, 'Value', 1);
    data.handles = handles;
    data = update_display(data, true);
  end

  % Release the GUI
  guidata(handles.hFig, data)
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function options_Callback(hObject, eventdata)
% This function is responsible for handling the buttons responsible for the
% option structure

  data = guidata(hObject);
  handles = data.handles;

  % Block the GUI
  set(handles.all_buttons, 'Enable', 'off');
  drawnow;
  refresh(hFig);

  % And get the type of button which called the callback (from its tag)
  type = get(hObject, 'tag');

  % By default, recompute
  recompute = true;

  % Handle all three buttons differently
  switch type

    % Save a snapshot
    case 'snapshot'

      % Fancy output
      disp('[Select a SVG filename]');

      % Prompting the user for the filename
      [fname, dirpath] = uiputfile({'*.svg', 'SVG vectorized image'}, ['Select a filename for your snapshot'], 'export/snapshot.svg');

      % Not cancelled
      if (ischar(fname))

        % This might take a while
        curr_name = get(handles.hFig, 'Name');
        set(hFig, 'Name', [curr_name ' (Saving snapshot...)']);

        % Get the full name and save the snapshot !
        fname = fullfile(dirpath, fname);
        plot2svg(fname, hFig);

        % And release !
        set(handles.hFig, 'Name', curr_name);
      end

      recompute = false;
  end

  % Release the GUI and recompute the display
  data = update_display(data, recompute);
  guidata(handles.hFig, data)
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function gui_Callback(hObject, eventdata)
% This function handles the callback of most buttons in the GUI !

  data = guidata(hObject);
  handles = data.handles;

  % By default we recompute the display
  recompute = true;

  % Get the channel index
  indx = handles.current;

  % If no more data, do nothing
  if (indx < 1)
    return;
  end

  % And get the type of button which called the callback (from its tag)
  type = get(hObject, 'tag');
  switch type

    % Each checkbox is responsible for its respective boolean fields
    case {'detrend', 'cosmics', 'hot_pixels', 'normalize'}
      data.channels(indx).(type) = logical(get(hObject, 'Value'));

    % A change in the channel index
    case 'channels'
      handles.current = get(hObject, 'Value');

    % Otherwise, do nothing. This is used to cancel the deletion requests
    otherwise
      return;
  end

  data.handles = handles;

  % Update the display accordingly
  data = update_display(data, recompute);

  guidata(handles.hFig, data)

  return
end

function cancel_CloseRequestFcn(hObject, eventdata)
% This function stops the current processing after confirmation

  data = guidata(hObject);
  handles = data.handles;

  % Just double check that the user want to quit
  answer = questdlg('Do you really want to discard all your changes ?');
  ok = strcmp(answer,'Yes');

  % If everything is OK, release the GUI and quit
  if (ok)
    delete(handles.hFig);
  end

  return
end

function channel_CloseRequestFcn(hObject, eventdata)
% This function converts the various indexes back into strings to prepare
% the channels structure for its standard form before releasing the GUI
% to exit it

  % If everything is OK, release the GUI and quit
  data = guidata(hObject);
  handles = data.handles;

  data = handle_rois(data, handles.current, 0);
  handles = data.handles;

  if (handles.current > 0)
    data.channels(handles.current).pixel_size = str2double(get(handles.resolution, 'String'));
    data.channels(handles.current).amplitude = str2double(get(handles.amplitude, 'String'));
  end

  % Copy the channels
  myrecording.channels = data.channels;
  % And get the experiment name
  myrecording.experiment = get(handles.experiment, 'String');

  [fpath, fname, fext] = fileparts(myrecording.channels(1).fname);
  mname = fullfile(fpath, [myrecording.experiment '.mat']);
  save(mname, 'myrecording');

  [junk, sizes] = find_zooids(myrecording);

  info = imfinfo(myrecording.channels(1).fname);

  if (exist('census.csv', 'file'))
    fid = fopen('census.csv', 'a');
  else
    fid = fopen('census.csv', 'a');
    fprintf(fid, 'Date, Name, Count\n');
  end
  fprintf(fid, '%s, %s, %d\n', info.FileModDate, mname, sum(sizes));
  fclose(fid);

  msgbox(['The segmentation and the analysis were saved to:\n' mname '\ncensus.csv']);

  delete(handles.hFig);

  return
end

function fnames = load_images()

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
  disp('[Select image(s) of the colony]');

  % Prompting the user for the movie file
  [fnames, dirpath] = uigetfile({'*.*'}, ['Load image(s) of the colony'], 'MultiSelect', 'on');
  fnames = fullfile(dirpath, fnames);

  % Return back to our original folder
  if(~isempty(curdir))
    cd(curdir);
  end

  % If no file was selected, stop here
  if (isempty(fnames)  ||  isequal(dirpath, 0))
    disp(['No image selected']);
    fnames = '';
    return;
  end

  %% All to cells
  if (ischar(fnames))
    fnames = {fnames};
  end

  return;
end
