function stk = filter_stack(file, alpha, Nsampling, border_erosion)

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

%%%