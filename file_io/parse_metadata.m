function [metadata, opts] = parse_metadata(data, opts)

  if (nargin < 2)
    opts = get_struct('options');
  end

  max_iter = 15;

  data = umanager2xml(data, max_iter);

  xml_data = parse_xml(data);

  xml_type = get_attribute(xml_data, 'xmlns');
  xml_type = regexp(xml_type, '^(http://)?(.*?)$', 'tokens');
  xml_type = xml_type{1}{2};

  switch xml_type
    case 'www.w3.org/uManager'
      summary_keys = {'Summary', 'Frames', 'Summary', 'Channels', 'Summary', 'Slices'};
      frame_keys = {'FrameKey*', 'Frame', 'FrameKey*', 'Slice', 'FrameKey*', 'Channel', 'FrameKey*', 'ElapsedTime-ms', 'FrameKey*', 'Exposure-ms', 'FrameKey*', 'Z-um'};
      infer_keys = {''};
      resol_keys = {''};
    case 'www.openmicroscopy.org/Schemas/OME/2012-06'
      summary_keys = {'Pixels', 'SizeT', 'SizeC', 'SizeZ'};
      frame_keys = {''};
      infer_keys = {''};
      resol_keys = {''};
    case 'tempuri.org/UltraviewSchm.xsd'
      summary_keys = {'Property', 'T', 'C', 'Z'};
      frame_keys = {''};
      infer_keys = {'ChannelSetting', 'ExposureTime','ChannelSetting', 'ChannelName',  'AcquiredTime', {'StartTime', 'FinishTime', 'yyyy-mm-ddTHH:MM:SS.FFF'}, 'ZSetting', {'TopPosition', 'BottomPosition'}};
      resol_keys = {'SpatialCalibration', {'XBasePixelSize', 'YBasePixelSize'}, 'SpatialCalibration', 'ObjectiveMagn', 'CameraSetting', 'Binning'};
    case 'schemas.datacontract.org/2004/07/LeicaMicrosystems.DataEntities.V3_2'
      summary_keys = {''};
      frame_keys = {'LasImage', '', '', '', 'AcquiredDate', '', 'Microscope_Focus_Position'};
      infer_keys = {'Camera', 'Exposure', 'Camera', 'Name', '', '', '', ''};
      resol_keys = {'LasImage', {'XMetersPerPixel', 'YMetersPerPixel'}, 'Microscope_Visual_Magnification', '', '', ''};
      keyboard
    otherwise
      summary_keys = {''};
      frame_keys = {''};
      infer_keys = {''};
      resol_keys = {''};

      warning(['Unknown XML schema "' xml_type '", unable to parse it.']);
  end

  keyboard

  metadata = parse_summary(xml_data, summary_keys);
  metadata = parse_frames(xml_data, frame_keys, metadata);
  metadata = infer_frames(xml_data, infer_keys, metadata);
  metadata.raw_data = data;

  opts = infer_resolution(xml_data, resol_keys, opts);

  return;
end

function value = get_attribute(node, attribute)

  value = '';
  if (isempty(node))
    return;
  end

  for i=1:length(node.Attributes)
    if (regexp(node.Attributes(i).Name, attribute))
      value = node.Attributes(i).Value;

      break;
    end
  end

  return;
end

function child = get_child(node, children)

  child = '';

  if (regexp(node.Name, children))
    child = node;
    return;
  end

  for i=1:length(node.Children)
    if (regexp(node.Children(i).Name, children))
      child = node.Children(i);

      break;
    end
  end

  if (isempty(child))
    for i=1:length(node.Children)
      child = get_child(node.Children(i), children);

      if (~isempty(child))

        break;
      end
    end
  end

  return;
end

function values = get_values(node, keys)

  nkeys = length(keys) / 2;

  values = cell(1, nkeys);

  for i=1:nkeys
    key = keys{2*i - 1};
    attr = keys{2*i};

    if (~isempty(key))
      params = {};
      if (iscell(key))
        params = key(2:end);
        key = key{1};
      end

      child = get_child(node, key);

      if (~isempty(child))
        if (isempty(attr))
          values{i} = child.Data;
        elseif (iscell(attr))
          nattrs = length(attr);
          tmp_val = cell(1, nattrs);

          for j=1:nattrs
            tmp_val{j} = get_attribute(child, attr{j});
          end
        else
          values{i} = get_attribute(child, attr);
        end
      end

      if (~isempty(params))
        values{i} = [values{i} params];
      end
    end
  end

  return;
end

function [channel_index, metadata] = get_channel(channel, metadata)

  is_empty = cellfun('isempty', metadata.channels);

  if (~all(is_empty))
    is_group = ismember(metadata.channels(~is_empty), channel);
  else
    is_group = ~is_empty;
  end

  if (any(is_group))
    channel_index = find(is_group, 1);
  else
    channel_index = find(is_empty, 1);

    metadata.channels{channel_index} = channel;
  end

  return;
end

function metadata = infer_frames(data, keys, metadata)

  nchar = length(keys{1});
  if (nchar == 0)
    return;
  end

  [nchannels, nframes, nslices] = size(metadata.acquisition_time);

  node = get_child(data, keys{1});
  for i=1:nchannels
    metadata.exposure_time(i, :, :) = str2double(get_attribute(node.Children(i).Children, keys{2}));
  end

  node = get_child(data, keys{3});
  for i=1:nchannels
    metadata.channels{i} = get_attribute(node.Children(i), keys{4});
  end

  node = get_child(data, keys{5});
  start = datevec(get_attribute(node, keys{6}{1}), keys{6}{3});
  ends = datevec(get_attribute(node, keys{6}{2}), keys{6}{3});
  step = 1000 * etime(ends, start) / (nframes + 1);

  for i=1:nframes
    metadata.acquisition_time(:, i, :) = (i-1)*step;
  end

  node = get_child(data, keys{7});
  top = str2double(get_attribute(node, keys{8}{1}));
  bottom = str2double(get_attribute(node, keys{8}{2}));
  step = (top - bottom) / max(nslices-1, 1);

  for i=1:nslices
    metadata.z_position(:, :, i) = top + (i-1)*step;
  end

  return;
end

function opts = infer_resolution(data, keys, opts)

  nchar = length(keys{1});
  if (nchar == 0)
    return;
  end

  node = get_child(data, keys{1});
  if (iscell(keys{2}))
    nres = length(keys{2});
    pixel_size = NaN(1, nres);
    for i=1:nres
      pixel_size(i) = str2double(get_attribute(node, keys{2}{i}));
    end
  else
    pixel_size = str2double(get_attribute(node, keys{2}));
  end

  node = get_child(data, keys{3});
  magnification = str2double(get_attribute(node, keys{4}));

  node = get_child(data, keys{5});
  binning = str2double(get_attribute(node, keys{6}));

  if (binning ~= 0 & magnification ~= 0 & pixel_size ~= 0)
    pixel_size = mean(pixel_size);

    opts.pixel_size = pixel_size;
    opts.magnification = magnification;
    opts.binning = binning;

    opts.ccd_pixel_size = pixel_size * magnification / binning;
  end

  return;
end


function metadata = parse_frames(data, keys, metadata)

  nchar = length(keys{1});
  if (nchar == 0)
    return;
  end

  nchild = length(data.Children);

  for i=1:nchild
    children = data.Children(i);

    values = get_values(children, keys);

    frame = str2double(values{1});
    slice = str2double(values{2});
    channel = values{3};

    if (isnan(frame))
      frame = 1;
    else
      frame = frame + 1;
    end
    if (isnan(slice))
      slice = 1;
    else
      slice = slice + 1;
    end

    [channel, metadata] = get_channel(channel, metadata);

    time = str2double(values{4});
    exposure = str2double(values{5});
    z_pos = str2double(values{6});

    metadata.acquisition_time(channel, frame, slice) = time;
    metadata.exposure_time(channel, frame, slice) = exposure;
    metadata.z_position(channel, frame, slice) = z_pos;
  end

  [junk, indexes] = sort(metadata.acquisition_time(:,1,1));

  metadata.channels = metadata.channels(indexes);
  metadata.acquisition_time = metadata.acquisition_time(indexes,:,:);
  metadata.exposure_time = metadata.exposure_time(indexes,:,:);
  metadata.z_position = metadata.z_position(indexes,:,:);

  return;
end

function metadata = parse_summary(xml_data, keys)

  metadata = get_struct('metadata');

  nchar = length(keys{1});
  if (nchar == 0)
    return;
  end

  values = get_values(xml_data, keys);

  nframes = str2double(values{1});
  nchannels = str2double(values{2});
  nslices = str2double(values{3});

  nframes = max(nframes, 1);
  nchannels = max(nchannels, 1);
  nslices = max(nslices, 1);

  data = NaN([nchannels, nframes, nslices]);

  metadata.channels = cell(nchannels, 1);

  metadata.acquisition_time = data;
  metadata.exposure_time = data;
  metadata.z_position = data;

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
