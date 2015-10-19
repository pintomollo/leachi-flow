function new_names = clean_borders(files)

  if nargin<1, files = '*.tif'; end

  new_names = {};
  dir_out = '_cleaned';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  fprintf(' Cleaning image borders :     ');

  hdil = strel('disk', 5);
  new_names = files;

  for i = 1:N % loop over images to resize images
    fprintf('\b\b\b%3d', i);

    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    orig = imread(filename);
    im = orig;

    c = size(im,3);

    borders = (any(isnan(im), 3) | all(im==0, 3));
    borders = ~(imfill(~borders, 'holes'));

    if (c > 1)
      im = rgb2gray(im);
    end

    im = adjust_tif(im);
    th = graythresh(im(:))*max(im(:));
    mask = (im > th);
    mask = imdilate(mask, hdil);
    mask = ~imfill(~(mask | borders), 'holes');
    mask = mask & ~borders;

    if (any(mask(:)))
      [im, bkgs] = imfillborder(orig);

      for j = 1:c
        tmp = im(:,:,j);
        tmp(mask) = bkgs(j);
        im(:,:,j) = tmp;
      end
      imwrite(im, new_name, 'TIFF');
    else
      copyfile(filename, new_name);
    end

    new_names{i} = new_name;
  end
  fprintf('\b\b\b\bdone\n');

  return;
end
