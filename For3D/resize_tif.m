function new_names = resize_tif(files, Amax, Bmax)

%% function resize_tif(files, Amax, Bmax)
%
% Put images at the same size, completing by black borders.
% Intended for building a stack from serial images in ImageJ.
%
% -- inputs --
% * files: generic name for all images (default: '*.tif')
% * Amax, Bmax: max number of lines & columns
% if not provided, or set to 0, max size is computed from the images
%
% -- output --
% Resized images are automatically saved in a subfolder named 'resized'.
%
% (c) AS 2013
%
% see also rendering_3D (For3D)

  if nargin<1, files = '*.tif'; end
  if nargin<3, Amax = 0; Bmax = 0; end

  new_names = {};
  dir_out = '_resized'; % '../resized';

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
      files{i} = fullfile(filepath, ls(i).name);
    end
  end

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  if (Amax == 0)
      A = zeros(N, 1); B = A;

      for i = 1:N % loop over images to find max size
          filename = files{i};
          im = imread(filename);
          [A(i), B(i), ~] = size(im); % ! lines & columns, not width & height !
          fprintf('%s, %i, %i\n', filename, A(i), B(i))
      end
  % % %     filenames = sort_nat(filenames);
      Amax = max(A);
      Bmax = max(B);
      fprintf('max size: %i, %i\n\n', Amax, Bmax)
  end

  new_names = files;

  for i = 1:N % loop over images to resize images
      filename = files{i};
      im = imread(filename);

      [a, b, c] = size(im);
      [filepath, fname, fileext] = fileparts(filename);
      out_path = fullfile(filepath, dir_out);

      if ~isdir(out_path)
        mkdir(out_path);
      end
      new_name = fullfile(out_path, [fname fileext]);

      fprintf('%s\n', new_name);

      if (~exist(new_name, 'file'))
        if (a < Amax || b < Bmax)

          a0 = round(Amax/2-a/2);
          b0 = round(Bmax/2-b/2);

          im2 = zeros(Amax, Bmax, c, class(im));
          im2(a0+1:a0+a, b0+1:b0+b, :) = im;
          im2 = imfillborder(im2);

          imwrite(im2, new_name, 'Compression', 'none')

        else
          copyfile(filename, new_name);
        end
      end
      new_names{i} = new_name;
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
