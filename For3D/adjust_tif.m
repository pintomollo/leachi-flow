function new_names = adjust_tif(files)

  if nargin<1, files = '*.tif'; end

  new_names = {};
  dir_out = '_adjusted'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Adjust images pixel values...');

  new_names = files;

  for i = 1:N % loop over images to resize images
      filename = files{i};
      im = imread(filename);

      [filepath, fname, fileext] = fileparts(filename);
      out_path = fullfile(filepath, dir_out);

      if ~isdir(out_path)
        mkdir(out_path);
      end
      new_name = fullfile(out_path, [fname fileext]);

      fprintf('%s\n', new_name);

      [im, bkgs] = imfillborder(im);

      %border_rows = [im(1,:,:) im(end,:,:)];
      %border_cols = [im(:,1,:); im(:,end,:)];

      is_white = (mean(bkgs) > mean(im(:)));
      if is_white
        mval = max(im(:));
        im = mval - im;
        bkgs = mval - bkgs;
      end % invert if white bkground

      for c = 1:size(im,3)
        tmp_img = im(:,:,c) - bkgs(c);
        tmp_img(tmp_img < 0) = 0;
        tmp_img = imnorm(tmp_img);
        im(:,:,c) = tmp_img;
      end

      imwrite(im, new_name, 'Compression', 'none')

      new_names{i} = new_name;
  end

  return;
end
