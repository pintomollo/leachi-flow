function params = filter_stack(params)

%% function stk = filter_stack(file, alpha, Nsampling, border_erosion)
%
% Filter image stack by equalizing fluo. levels, 2D median and 3D gaussian filters.
%
% -- inputs --
% * file: image stack filename
% * alpha: coefficient for the sub-function equalize_stack (from 0 to 1)
% * Nsampling: number for down-sampling images
% * border_erosion: number of pixels to be removed (intensity set to 0)
% in the red & green channels, respectively to the thresholded blue channel
%
% -- output --
% Filtered images are automatically saved, with the extension '_filt.tif'.
%
% (c) Arnauld SERGE, 2011
%
% see also rendering_3D, equalize_stack

  if nargin < 1
    params = get_struct('For3D');
  end

  files = params.filename;

  stk_files = {};
  dir_out = '_tmp';

  [files, out_path] = get_filenames(files, dir_out);

  Nz = length(files);
  if (Nz == 0), disp('nada??'), return, end
  %Nk = length(imfinfo(file));
  %kk = (1:Nk);

  new_names = files;

  marg = 3;
  %Msize = [7 7]; % odd size is better
  Msize = 3;
  Gsd = 0.65; % Gsize = 5; Gsd = 2

  borders = any(params.border_erosion);

  if (borders)
    Nc = length(params.border_erosion);
    disk = cell(Nc,1);
    for nc = 1:Nc
        disk{nc} = strel('disk', round(border_erosion(nc)));
    end
  end

  fprintf(' Median filtering over %i pixels radius & gaussian blur of %f sigma for %i images :     ', Msize, Gsd, Nz)
  for nz = 1:Nz
    fprintf('\b\b\b%3d', nz);

    filename = files{nz};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    im = imread(filename);
    type = class(im);

    im = double(im);

    imsize = size(im);
    imsize(1:2) = imsize(1:2) + 2*marg;
    tmp_img = zeros(imsize);
    tmp_img(1+marg:end-marg, 1+marg:end-marg, :) = im;
    im = tmp_img;

    clear tmp_img;

    im = median_mex(im, Msize);
    im = gaussian_mex(im, Gsd);

    if (borders)
      if (size(im, 3) > 1)
        ref = rgb2gray(im);
      else
        ref = img;
      end

      if (max(ref(:)) > 0)
        ref = imfill(ref);
        ref = (ref > mean(ref(:)));

        for nc = 1:size(im, 3)
          if nc <= length(border_erosion) && border_erosion(nc) > 0
            tmp_ref = imerode(ref, disk{nc});
            im_c = im(:, :, nc);
            im_c(tmp_ref==0) = 0;
            im(:,:,nc) = im_c;
          end
        end

        clear tmp_ref im_c;
      end

      clear ref;
    end

    imwrite(cast(im, type), new_name, 'TIFF');

    new_names{nz} = new_name;
  end

  clear im;

  fprintf('\b\b\b\bdone\n');

  tmp_stack = write_stack(new_names);
  dir_out = '_filtered';

  tmp_dir = out_path;
  [tmp_stack, out_path] = get_filenames(tmp_stack, dir_out);

  Nc = length(tmp_stack);
  threshs = zeros([1 Nc]);

  if (isempty(params.sparse_thresholds))
    params.sparse_thresholds = NaN([1, Nc]);
  elseif (numel(params.sparse_thresholds) ~= Nc)
    tmp = params.sparse_thresholds(1:min(Nc,end));
    tmp = tmp(:).';
    tmp = [tmp NaN([1, Nc-length(tmp)])];
    params.sparse_thresholds = tmp;
  end

  fprintf(' Z-gaussian blur of %f sigma for %i channels and size min %f :     ', Gsd, Nc, params.min_volume);
  for i=1:Nc
    fprintf('\n  %d :', i);

    filename = tmp_stack{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    [stk, threshs(i), type] = load_sparse_stack(filename, params.sparse_thresholds(i));
    stk = imssmooth(stk, 3, Gsd, threshs(i)/2);

    cc = bwconncomp(full(stk>0));
    min_volume = floor(prod(cc.ImageSize) / (params.min_volume^3));

    for j=1:cc.NumObjects
      pix = cc.PixelIdxList{j};
      if (length(pix) <= min_volume)
        stk(pix) = 0;
      end
    end

    write_stack(new_name, stk, type, num2str(threshs(i)));
    stk_files{i} = new_name;

  end

  [split_dir, junk, junk] = fileparts(tmp_stack{1});

  rmdir(tmp_dir, 's');
  rmdir(split_dir, 's');

  fprintf('\n  done !\n');

  params.sparse_thresholds = threshs;
  params.filename = stk_files;
%%%
