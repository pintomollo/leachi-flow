function resize_tif(files, Amax, Bmax)

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

dir_out = 'resized'; % '../resized';
if ~isdir(dir_out), mkdir(dir_out), end

files = dir(files);
N = length(files);
if (N == 0), disp('nada??'), return, end

if (Amax == 0)
    A = zeros(N, 1); B = A;
    
    for i = 1:N % loop over images to find max size
        filename = files(i).name;
        im = imread(filename);
        [A(i), B(i), ~] = size(im); % ! lines & columns, not width & height !
        fprintf('%s, %i, %i\n', filename, A(i), B(i))
    end
% % %     filenames = sort_nat(filenames);
    Amax = max(A);
    Bmax = max(B);
    fprintf('max size: %i, %i\n\n', Amax, Bmax)
end

if ~isdir(dir_out), mkdir(dir_out), end

% b_max = 3000 % largest allowed? cf variable zoom in calopix

for i = 1:N % loop over images to resize images
    filename = files(i).name;
    im = imread(filename);
    
    [a, b, c] = size(im);
%     if b>b_max, im = imresize(im, 1/2); [a, b, c] = size(im); fprintf('image too large, %i pxl width, resized\r', b), end
    
    border_rows = squeeze([im(1,:,:) im(end,:,:)]);
    border_cols = squeeze([im(:,1,:); im(:,end,:)]);
    background_val = mean([border_rows; border_cols]);
    
    fprintf('%s\n', filename)
    a0 = round(Amax/2-a/2);
    b0 = round(Bmax/2-b/2);
    
    im2 = uint8(zeros(Amax, Bmax, c));
    for j = 1:c
        im2(:,:,j) = im2(:,:,j) + background_val(j);
    end
    im2(a0+1:a0+a, b0+1:b0+b, :) = im;
    imwrite(im2, [dir_out filesep filename], 'Compression', 'none')
end

%%%