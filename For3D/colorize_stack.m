function new_names = colorize_stack(files)

  if nargin<1, files = '*.tif'; end

  new_names = {};
  dir_out = '_colored'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Coloring images...');

  new_names = files;

  filename = files{1};
  im = imread(filename);
  [junk, nhist] = imsplitcolors(im);

  for i = 2:N % loop over images to resize images
    filename = files{i};
    im = imread(filename);

    [junk, nhist] = imsplitcolors(im, nhist);
  end

  [imax] = imsplitcolors(nhist);

  for i = 1:N % loop over images to resize images
      filename = files{i};
      im = imread(filename);

      [filepath, fname, fileext] = fileparts(filename);
      new_name = fullfile(out_path, [fname fileext]);

      fprintf('%s\n', new_name);

      im = imsplitcolors(im, imax);

      imwrite(im, new_name, 'Compression', 'none');

      new_names{i} = new_name;
  end

  return;
end
