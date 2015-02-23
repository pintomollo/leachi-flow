function [metadata, opts] = parse_metadata(data, opts)

  if (nargin < 2)
    opts = get_struct('options');
  end

  metadata = get_struct('metadata');
  max_iter = 15;

  data = umanager2xml(data, max_iter);
  data = parse_xml(data);

  keyboard

  return;

  data = strsplit(data, {'\n', '\r'});
  data = data(~cellfun('isempty', data));

  if (~isempty(data))

    size_increment = 200;
    vector_size = size_increment;

    nans = NaN(1, size_increment);
    pos = nans;
    frame = nans;
    plane = nans;
    group = nans;
    exposure = nans;
    time = nans;
    groups = {};
    pixel_size = 0;
    binning = 0;
    magnification = 0;

    count = 1;
    state = 'init';
    block_level = 0;
    obj_level = 0;

    said_it = false;

    % Read the first line
    line = data{1};

    if (any(line == '{'))
      % uManager metadata file

      % Loop until we have read the whole file
      for i=1:length(data)
        line = data{i};

        if (any(line == '{'))
          block_level = block_level + 1;
        end
        if (any(line == '}'))
          block_level = block_level - 1;
        end

        switch state
          case 'init'
            if (any(line == '{'))
              tokens = regexp(line, '^\s*"(.+)": {$','tokens');
              if (isempty(tokens))
                state = 'init';
              elseif (strncmp(tokens{1}, 'FrameKey', 8))
                tmp_tokens = regexp(tokens{1}, 'FrameKey-(\d+)-\d+-\d+','tokens');
                if (~isempty(tmp_tokens) & str2double(tmp_tokens{1}{1}) >= nframes)
                  if (~said_it)
                    warning(['Metadata describes a plane (' tmp_tokens{1}{1}{1} ') that is out of the range of the actual recordings (' num2str(nframes-1) '), ignoring this infromation.']);
                    said_it = true;
                  end

                  state = 'ignore';
                  obj_level = block_level;
                else
                  state = 'read';
                  obj_level = block_level;
                end
              else
                state = 'ignore';
                obj_level = block_level;
              end
            end
          case 'ignore'
            if (~isempty(strfind(line, 'PixelSize')))
              tokens = regexp(line, '^\s*"(.+)": "?(.*?)"?,?$','tokens');
              pixel_size = str2double(tokens{1}{2});
            end
            if (block_level < obj_level)
              state = 'init';
            end
          case 'read'
            if (block_level < obj_level)
              count = count + 1;
              if (count > vector_size)
                vector_size = vector_size + size_increment;
                pos = [pos nans];
                frame = [frame nans];
                plane = [plane nans];
                group = [group nans];
                exposure = [exposure nans];
                time = [time nans];
              end

              state = 'init';
            else
              tokens = regexp(line, '^\s*"(.+)": "?(.*?)"?,?$','tokens');
              if (~isempty(tokens))
                tokens = tokens{1};
              end
              
              if (numel(tokens) == 2)
                switch tokens{1}
                  case 'Z-um'
                    pos(count) = str2double(tokens{2});
                  case 'Frame'
                    frame(count) = str2double(tokens{2}) + 1;
                  case 'Slice'
                    plane(count) = str2double(tokens{2}) + 1;
                  case 'Channel'
                    is_group = ismember(groups, tokens{2});
                    if (isempty(is_group)|~any(is_group))
                      groups{end+1} = tokens{2};
                      is_group = [is_group true];
                    end
                    group(count) = find(is_group);
                  case 'Exposure-ms'
                    exposure(count) = str2double(tokens{2});
                  case 'ElapsedTime-ms'
                    time(count) = str2double(tokens{2});
                  case {'Time', 'FileName'}
                    % Ignoring
                  otherwise
                    disp(['Parsing incomplete ...']);
                end
              end
            end
          otherwise
            disp(['Unkown parsing state ("' state '"), skipping...']);
            state = 'init'; 
        end
      end

      goods = (~isnan(frame) & ~isnan(plane));
      frame = frame(goods);
      plane = plane(goods);
      pos = pos(goods);
      time = time(goods);
      exposure = exposure(goods);
      group = group(goods);

      frames = unique(frame);
      planes = unique(plane);
      channels = unique(group);
      sizes = [length(channels), length(frames), length(planes)];

      order = sub2ind(sizes(2:3), frame, plane);
      sizes = [sizes(1) prod(sizes(2:3))];
      order = sub2ind(sizes, group, order);

      tmp = NaN(sizes);
      tmp(order) = pos;
      pos = tmp;
      tmp(order) = frame;
      frame = tmp;
      tmp(order) = plane;
      plane = tmp;
      tmp(order) = time;
      time = tmp;
      tmp(order) = exposure;
      exposure = tmp;
      tmp(order) = group;
      group = tmp;

      [iindx, jindx] = find(isnan(pos));

      for i = iindx
        for j = jindx
          if (j > 1)
            pos(i, j) = pos(i, j-1);
            frame(i, j) = frame(i, j-1) + 1;
            plane(i, j) = plane(i, j-1);
            time(i, j) = time(i, j-1);
            exposure(i, j) = exposure(i, j-1);
            group(i, j) = group(i, j-1);
          end
        end
      end

    elseif (any(line == '<'))
      % Spinning Disk metadata file

      % Loop until we have read the whole file
      for i=1:length(data)
        line = data{i};

        if (~isempty(strfind(line, '</')))
          if (sum(line == '<') == 1)
            block_level = block_level - 1;
          end
        elseif (isempty(strfind(line, '/>')))
          block_level = block_level + 1;
        end

        switch state
          case 'init'
            if (~any(line == '/'))
              tokens = regexp(line, '^\s*<(.+)>$','tokens');
              if (isempty(tokens))
                state = 'init';
              elseif (strncmp(tokens{1}, 'CameraSetting', 13)) | (strncmp(tokens{1}, 'SpatialCalibration', 18))
                state = 'read';
                obj_level = block_level;
              end
            end
          case 'read'
            if (block_level < obj_level)

              state = 'init';
            else
              tokens = regexp(line, '^\s*<(.+)>(.*?)</\1>$','tokens');
              if (~isempty(tokens))
                tokens = tokens{1};
              end
              
              if (numel(tokens) == 2)
                switch tokens{1}
                  case 'ObjectiveMagn'
                    magnification = str2double(tokens{2});
                  case 'XBasePixelSize'
                    pixel_size(1) = str2double(tokens{2});
                  case 'YBasePixelSize'
                    pixel_size(2) = str2double(tokens{2});
                  case 'Binning'
                    binning = str2double(tokens{2});
                  case 'ExposureTime'
                    exposure = str2double(tokens{2});
                  otherwise
                    %disp(['Parsing incomplete ...']);
                end
              end
            end
          otherwise
            disp(['Unkown parsing state ("' state '"), skipping...']);
            state = 'init'; 
        end
      end
    else
      warning('Unknown metadata file type, ignoring it.');

      return;
    end

    metadata.z_position = pos;
    metadata.frame_index = frame;
    metadata.plane_index = plane;
    metadata.acquisition_time = time;
    metadata.exposure_time = exposure;
    metadata.channel_index = group;
    metadata.channels = groups;

    if (binning ~= 0 & magnification ~= 0 & pixel_size ~= 0)
      pixel_size = mean(pixel_size);

      opts.pixel_size = pixel_size;
      opts.magnification = magnification;
      opts.binning = binning;

      opts.ccd_pixel_size = pixel_size * magnification / binning;
    end
  end

  return;
end

function data = umanager2xml(data, max_iter)

  if (any(data == '{'))
    rem_spaces = @remove_spaces;
    rem_eol = @remove_eol;

    data = regexprep(data, '("\S*? +\S*?":)', '${rem_spaces($1)}');
    data = regexprep(data, '(\[[^\]]*?\],)', '${rem_eol($1)}');

    data = regexprep(data, '"([^\n\r]*?)": ([^\n\r{]*?),?[\r\n]+', '<$1>$2</$1>\n');
    nprev = length(data);
    for i=1:max_iter
      data = regexprep(data, '[\n\r](\s+)"([^\n\r]*?)": {(.*?)[\n\r]\1}[,\r\n]+', '\n$1<$2>$3\n$1</$2>\n');
      curr = length(data);
      if (nprev == curr)
        break;
      else
        nprev = curr;
      end
    end
    data = regexprep(data, '{(.*)}', ['<?xml version="1.0" encoding="utf-8"?>\n'...
                                        '<TimeSeries xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '...
                                        'xmlns="http://www.w3.org/uManager">$1</TimeSeries>']);
    data = regexprep(data, '>"', '>');
    data = regexprep(data, '"<', '<');
  end

  return;
end

function str = remove_eol(str)

  str = regexprep(str, '\s*[\n\r]\s*', '');

  return;
end

function str = remove_spaces(str)

  str = strrep(str, ' ', '_');

  return;
end
