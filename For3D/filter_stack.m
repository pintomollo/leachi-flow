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

  fprintf('\n Median filtering over %i pixels rdius & gaussian blur of %f sigma for %i images :    ', Msize, Gsd, Nz)
  for nz = 1:Nz
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
    fprintf('\b\b\b%3d', nz)
  end

  clear im;

  tmp_stack = write_stack(new_names);
  dir_out = '_filtered';

  tmp_dir = out_path;
  [tmp_stack, out_path] = get_filenames(tmp_stack, dir_out);

  Nc = length(tmp_stack);
  threshs = zeros([1 Nc]);

  if (isempty(params.sparse_thresholds))
    params.sparse_thresholds = NaN([1, Nc]);
  end

  fprintf('\n Z-gaussian blur of %f sigma for %i channels :    ', Gsd, Nc);
  for i=1:Nc
    filename = tmp_stack{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    [stk, threshs(i), type] = load_sparse_stack(filename, params.sparse_thresholds(i));
    stk = imssmooth(stk, 3, Gsd, threshs(i)/2);

    write_stack(new_name, stk, type, num2str(threshs(i)));
    stk_files{i} = new_name;

    fprintf('\b\b\b%3d', i)
  end

  rmdir(tmp_dir, 's');

  params.sparse_thresholds = threshs;
  params.filename = stk_files;

%{
smooth_file = ['temp' filesep file(1:end-4) '_smoothed.tif'];
im = imread(smooth_file, 'Index', 1);
Nx = size(im, 1);
Ny = size(im, 2);
Nz = length(imfinfo(smooth_file));
Nc = size(im, 3);

marg = 3;
Nx_out = ceil(Nx/Nsampling) + 2*marg;
Ny_out = ceil(Ny/Nsampling) + 2*marg;
Nz_out = Nz + 2*marg;
stk = zeros(Nx_out, Ny_out, Nz_out, Nc);
eq_file = ['temp' filesep file(1:end-4) '_eq.tif'];
filt_file = [file(1:end-4) '_filt.tif'];


if isempty(dir(eq_file))
    
    %% ****** data: load, equalize, smooth (median filter size 8 pxl) & redim (subsampling in x & y) ******
    fprintf(' Starting %s, with %d images \n', smooth_file, Nz)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if alpha>0, equalize_stack(smooth_file, alpha), end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if isempty(dir(filt_file))
    
    Msize = [7 7]; % odd size is better
    
    fprintf('\n Median filtering over %ix%i pixels & sampling by a factor of %i %s, image    ', Msize, Nsampling, file)
    for nz = 1:Nz
        if alpha>0, im_z = imread(eq_file, 'Index', nz);
        else im_z = imread(smooth_file, 'Index', nz);
        end
        
        %% median filter 7x7 (images, 2D) & sampling
        for nc = 1:Nc
            im_filt = medfilt2(im_z(:, :, nc), Msize);
            stk(1+marg:end-marg, 1+marg:end-marg, nz+marg, nc) = imresize(im_filt, 1/Nsampling);
        end
        fprintf('\b\b\b%3d', nz)
    end
    
    %% filter by Gaussian (volume, 3D)
    Gsize = 3; Gsd = 0.65; % Gsize = 5; Gsd = 2
    
    fprintf('\n Gaussian filtering in 3D, over %ix%ix%i pixels, with %g sd...  ', Gsize, Gsize, Gsize, Gsd)
    for nc = 1:Nc
        stk(:, :, :, nc) = smooth3(stk(:, :, :, nc), 'gaussian', Gsize, Gsd);
        % % %         %% background subtraction for local (adaptative) thresholding
        % % %         GsizeB = 31; GsdB = 15;
        % % %         for nz = 1:Nz
        % % %             nnz = max(nz-1,1):min(nz+1,Nz);
        % % %             Bkgnd =  smooth3(stk(:, :, nnz, nc), 'gaussian', GsizeB, GsdB);
        % % %             stk(:, :, nnz, nc) = stk(:, :, nnz, nc) - Bkgnd + mean(Bkgnd(:));
        % % %             fprintf('\b\b\b%3d', nz)
        % % %         end
    end
    
    disk = cell(Nc,1);
    for nc = 1:Nc-1
        disk{nc} = strel('disk', round(border_erosion(nc)/Nsampling));
    end
    
    fprintf('\n Saving filtered image    ')
    for nz = 1:Nz
        im_filt = stk(1+marg:end-marg, 1+marg:end-marg, nz+marg, :);
        im_filt = squeeze(im_filt);
        
        %% remove border, ie capsule for possible artefactual labelling in red
        Blue = im_filt(:, :, 3); % blue image at nz, ref for erosion
        if (Nc > 1) && (max(Blue(:)) > 0)
            Blue = imfill(Blue);
            Blue = Blue > mean(Blue(:));
            for nc = 1:Nc-1
                if border_erosion(nc) > 0
                    Blue = imerode(Blue, disk{nc});
                    im_c = im_filt(:, :, nc);
                    im_c(Blue==0) = 0;
                    im_filt(:, :, nc) = im_c;
                end
            end
        end
        
        if (nz==1), imwrite(uint8(im_filt), filt_file, 'tiff', 'Compression', 'none')
        else imwrite(uint8(im_filt), filt_file, 'tiff', 'Compression', 'none', 'WriteMode', 'append')
        end
        fprintf('\b\b\b%3d', nz)
    end
    
else %% loading already filtered images
    fprintf('loading filtered image for %s    ', filt_file)
    for nz = 1:Nz
        im_z = imread(filt_file, 'Index', nz);
        if Nc == 1, im_z = im_z(:, :, 1); end
        if size(stk(1+marg:end-marg, 1+marg:end-marg, nz+marg, :), 1) ~= size(im_z, 1)
            disp('Wrong size for sampled image in altreasy existing filt file, please check pixel size or delete old filt file')
            stk = []; return
        end
        stk(1+marg:end-marg, 1+marg:end-marg, nz+marg, :) = im_z;
        fprintf('\b\b\b%3d', nz)
    end
end % filt

fprintf('\n')
%}
%%%
