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

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  %disp('Resizing images...');
  fprintf(' Parsing images to find the maximum size : ');

  pmsg = '';
  if (Amax == 0)
      A = zeros(N, 1); B = A;

      for i = 1:N % loop over images to find max size
          filename = files{i};
          im = imread(filename);
          [A(i), B(i), ~] = size(im); % ! lines & columns, not width & height !
          %fprintf('%s, %i, %i\n', filename, A(i), B(i))
          msg = sprintf('%i x %i', A(i), B(i));
          fprintf([repmat('\b', 1, length(pmsg)) repmat(' ', 1, length(pmsg)) repmat('\b', 1, length(pmsg))]);
          fprintf(msg);
          pmsg = msg;
      end
  % % %     filenames = sort_nat(filenames);
      Amax = max(A);
      Bmax = max(B);
      %fprintf('max size: %i, %i\n\n', Amax, Bmax)

      msg = sprintf('%i x %i is the max size\n', Amax, Bmax);
      fprintf([repmat('\b', 1, length(pmsg)) msg]);
  end

  new_names = files;

  fprintf(' Resizing images to the max size :     ');
  for i = 1:N % loop over images to resize images
      fprintf('\b\b\b%3d', i);

      filename = files{i};
      im = imread(filename);

      [a, b, c] = size(im);
      [filepath, fname, fileext] = fileparts(filename);

      new_name = fullfile(out_path, [fname fileext]);

      %fprintf('%s\n', new_name);

      if (~exist(new_name, 'file'))
        if (a < Amax || b < Bmax)

          a0 = round(Amax/2-a/2);
          b0 = round(Bmax/2-b/2);

          im2 = zeros(Amax, Bmax, c, class(im));
          im2(a0+1:a0+a, b0+1:b0+b, :) = im;

          [im2] = imfillborder(im2);
          imwrite(im2, new_name, 'Compression', 'none')

        else
          copyfile(filename, new_name);
        end
      end
      new_names{i} = new_name;
  end
  fprintf('\b\b\b\bdone\n');

  return;
end
