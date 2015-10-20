function [new_names, imax] = colorize_stack(files, colors)

  if nargin<1, files = '*.tif'; end
  if nargin<2, colors = []; end

  new_names = {};
  dir_out = '_colored'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  new_names = files;

  if (isempty(colors) || any(~isfinite(colors(:))) || ~isnumeric(colors))
    %disp('Coloring images...');
    fprintf(' Computing the color panel of the images :     ');

    new_names = files;


    ndist = 10;
    nhist = zeros(255);

    %filename = files{1};
    %im = imread(filename);
    %[junk, nhist] = imsplitcolors(im);

    for i = 1:N % loop over images to resize images
      fprintf('\b\b\b%3d', i);

      filename = files{i};
      im = imread(filename);

      [junk, nhist] = imsplitcolors(im, ndist, nhist);
    end

    [imax] = imsplitcolors(nhist, ndist);

    fprintf('\b\b\b\bdone\n');
  else
    imax = colors;
  end

  fprintf(' Coloring the images :     ');

  for i = 1:N % loop over images to resize images
      fprintf('\b\b\b%3d', i);

      filename = files{i};
      im = imread(filename);

      [filepath, fname, fileext] = fileparts(filename);
      new_name = fullfile(out_path, [fname fileext]);

      %fprintf('%s\n', new_name);

      im = imsplitcolors(im, imax);

      imwrite(im, new_name, 'Compression', 'none');

      new_names{i} = new_name;
  end

  fprintf('\b\b\b\bdone\n');

  return;
end
