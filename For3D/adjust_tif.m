function new_names = adjust_tif(files)

  if nargin<1, files = '*.tif'; end

  if (isnumeric(files))

    new_names = adjust_img(files);
  else

    new_names = {};
    dir_out = '_adjusted'; % '../resized';

    [files, out_path] = get_filenames(files, dir_out);

    N = length(files);
    if (N == 0), disp('nada??'), return, end

    fprintf(' Adjusting the pixel values of the images :     ');

    new_names = files;

    for i = 1:N % loop over images to resize images
        fprintf('\b\b\b%3d', i);

        filename = files{i};
        im = imread(filename);

        [filepath, fname, fileext] = fileparts(filename);
        new_name = fullfile(out_path, [fname fileext]);

        im = adjust_img(im);

        imwrite(im, new_name, 'Compression', 'none');

        new_names{i} = new_name;
    end

    fprintf('\b\b\b\bdone\n');
  end

  return;
end

function img = adjust_img(img)

  [img, bkgs] = imfillborder(img);

  %border_rows = [im(1,:,:) im(end,:,:)];
  %border_cols = [im(:,1,:); im(:,end,:)];

  is_white = (mean(bkgs) > mean(img(:)));
  if is_white
    mval = max(img(:));
    img = mval - img;
    bkgs = mval - bkgs;
  end % invert if white bkground

  for c = 1:size(img,3)
    tmp_img = img(:,:,c) - bkgs(c);
    tmp_img(tmp_img < 0) = 0;
    tmp_img = imnorm(tmp_img);
    img(:,:,c) = tmp_img;
  end

  return;
end
