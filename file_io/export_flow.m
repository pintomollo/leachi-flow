function export_flow(myrecording, nstep, opts)

  % Input checking and default values
  if (nargin < 1)
    return;
  elseif (nargin < 2)
    opts = get_struct('options');
    nstep = 10;
  elseif (nargin < 3)
    if (isstruct(nstep))
      opts = nstep;
      nstep = 10;
    else
      opts = get_struct('options');
    end
  end

  fname = myrecording.experiment;

  % Check if there is a folder name in the name itself
  [filepath, name, ext] = fileparts(fname);

  % Otherwise, put them into the export directory
  if (isempty(filepath))
    filepath = 'export';
  end

  % Create the directory
  if (~exist(filepath, 'dir'))
    mkdir(filepath);
  end

  % Build the full name
  fname = fullfile(filepath, name);


  % Obtain the pixel size of the screen
  set(0,'units','pixels');
  screen_size = get(0,'screensize');


  % And the figure name
  fig_name = 'Recording video, do not hide this window !! ';

  % Prepare the figure and the axis
  hFig = figure('Visible', 'off', ...
                'NumberTitle', 'off', ...
                'Units', 'pixels', ...
                'Name', fig_name);

  hAxes = axes('Parent', hFig, ...
               'DataAspectRatio', [1 1 1], ...
                'Units', 'pixels', ...
               'Visible', 'off',  ...
               'NextPlot', 'add', ...
               'Tag', 'axes');

  hPlot = axes('Parent', hFig, ...
               'Visible', 'off',  ...
                'Units', 'pixels', ...
               'NextPlot', 'add', ...
               'Tag', 'axes');

  % Set up the handlers
  hImg = -1;
  hBranch = -1;
  hFlow = -1;
  hNow = -1;

  % Initialize some computer-specific variables required for the conversion
  maxuint = intmax('uint16');

  % Open the specified AVI file with the maximal quality
  movie_name = [fname '.avi'];
  mymovie = VideoWriter(movie_name);
  mymovie.Quality = 100;
  open(mymovie);

  % Get the current number of frames and format the corresponding part of the title
  nframes = size_data(myrecording.channels(1));
  total_str = ['/' num2str(nframes)];

  vasculature = myrecording.segmentations.detections(1).properties;

  x1 = vasculature(1:3:end, 1);
  x2 = vasculature(2:3:end, 1);
  y1 = vasculature(1:3:end, 2);
  y2 = vasculature(2:3:end, 2);
  all_signs = vasculature(1:3:end, 4);

  valids = (all_signs ~= 0);
  x1 = x1(valids);
  x2 = x2(valids);
  y1 = y1(valids);
  y2 = y2(valids);
  all_signs = all_signs(valids);

  centers = 0.5*[(x1+x2) (y1+y2)];
  vects = [x2-x1 y2-y1];
  norms = bsxfun(@rdivide, vects, sqrt(sum(vects.^2, 2)) .* all_signs);

  flow_params = myrecording.trackings.detections.properties;

  timepts = [1:nframes]*opts.time_interval;
  speeds = opts.time_interval*cossum(timepts, flow_params)/opts.pixel_size;

  % Loop over all frames
  for n = 1:nstep:nframes

    nimg = timepts(n);
    speed = speeds(n)*nstep;

    % Get the image
    img = double(load_data(myrecording.channels(1), n));

    % We either replace the image or create a new one
    if (ishandle(hImg))
      set(hImg, 'CData', img);
    else

      % We need the image size to set the axes properly
      ssize = size(img);

      maxsize = 1 / max(ssize .* [1.6 1] ./ screen_size([4 3]));

      % Adapt the size of the image, fix the aspect ratio and the pixel range
      set(hFig, 'Visible', 'on', 'Position', [1 1 ssize([2 1]) .* [1 1.6] * maxsize], ...
          'MenuBar', 'none', 'ToolBar', 'none');
      hImg = image(img,'Parent', hAxes, 'CDataMapping', 'scaled');
      set(hAxes,'Visible', 'off', 'CLim', [0 maxuint], 'Position', [1 1 ssize([2 1]) * maxsize], ...
          'XLim', [1 ssize(2)], 'YLim', [1 ssize(1)], 'DataAspectRatio',  [1 1 1], ...
          'YDir', 'reverse');
      set(hPlot, 'XLim', [timepts([1 end])], 'YLim', [-1 1]*max(abs(speeds)));
      border = get(hPlot, 'TightInset') + 3;
      set(hPlot, 'Visible', 'on', 'Position', [border(1) (ssize(1) * maxsize)+border(2) (ssize(2) * maxsize)-border(1) (ssize(1) * maxsize)/2-border(2)]);


      % Use the defined colormap
      colormap(hFig, gray);

    end

    if (~ishandle(hBranch))
      set(hAxes,'NextPlot', 'add');
      hBranch = line('XData', vasculature(:,1), 'YData', vasculature(:,2), 'Parent', hAxes, 'Color', 'b', 'Marker', 'o');
    end

    if (ishandle(hFlow))
      set(hFlow, 'UData', norms(:,1)*speed, 'VData', norms(:,2)*speed);
    else
      hFlow = quiver(centers(:,1), centers(:,2), norms(:,1)*speed, norms(:,2)*speed, 0, 'Parent', hAxes, 'Color', 'r', 'LineWidth', 1.5);
    end

    if (ishandle(hNow))
      set(hNow, 'XData', [nimg nimg]);
    else
      line('XData', get(hPlot,'XLim'), 'YData', [0 0], 'Parent', hPlot, 'Color', 'k');
      line('XData', timepts, 'YData', speeds, 'Parent', hPlot, 'Color', 'r', 'LineWidth', 1.5);
      hNow = line('XData', [nimg nimg], 'YData', get(hPlot, 'YLim'), 'Parent', hPlot, 'Color', 'b');
    end

    % Update the name as a status bar
    set(hFig, 'Name', [fig_name num2str(n) total_str])

    % Get the current frame and store it in the movie
    frame = getframe(hFig);
    writeVideo(mymovie, frame);
    drawnow
  end

  % Close the movie
  close(mymovie)

  % Delete the figure
  delete(hFig)


  return
end
