function [files, out_path] = get_filenames(files, dir_out)

  out_path = '';

  if nargin<1
    files = '*.tif';
    dir_out = '';
  end

  if nargin < 2
    dir_out = '';
  end

  if (iscell(files) && length(files)==1)
    files = files{1};
  end

  if (~iscell(files))
    [filepath, filename, fileext] = fileparts(files);
    ls = dir(files);
    ls = clean_dir(ls);

    N = length(ls);

    if (N == 0)
      ls = dir(fullfile(files, '*.tif'));
      ls = clean_dir(ls);
      N = length(ls);
    end

    files = cell([N 1]);

    for i = 1:N
      files{i} = ls(i).name;
    end
    files = natsort(files);
    for i = 1:N
      files{i} = absolutepath(fullfile(filepath, files{i}));
    end
  end

  if (nargout > 1 && ~isempty(dir_out) && ~isempty(files))
    [filepath, fname, fileext] = fileparts(files{1});
    [shorter_path, prev_dir] = fileparts(filepath);

    if (prev_dir(1) == '_')
      out_path = fullfile(shorter_path, dir_out);
    else
      out_path = fullfile(filepath, dir_out);
    end

    if ~isdir(out_path)
      mkdir(out_path);
    end
  end

  return;
end

function list = clean_dir(list)

  list = list(~[list.isdir]);
  for i=length(list):-1:1
    if (list(i).name(1)=='.')
      list(i) = [];
    end
  end

  return;
end
