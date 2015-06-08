function fill_corners(file)

%% function fill_corners(file)
%
% Find black values and fill the with background color.
%
% -- input --
% file: image stack filename
%
% -- output --
% For an input stack named 'input_file.tif', filled images are automatically
% saved in a subfolder named 'input_file__black_corners_filled'.
%
% (c) Arnauld SERGE, 2011
%
% see also rendering_3D, equalize_stack


if nargin<1, files = dir('*.tif'); file = files(1).name; end

dir_out = [file '_black_corners_filled'];
if ~isdir(dir_out), mkdir(dir_out), end

N = length(imfinfo(file));

for i = 1:N
    im = imread(file, i);
    
    border_rows = squeeze([im(1,:,:) im(end,:,:)]); % border pixels
    border_cols = squeeze([im(:,1,:); im(:,end,:)]);
    background_val = mean([border_rows(sum(border_rows,2)>0,:); border_cols(sum(border_cols,2)>0,:)]); % mean color of border (non black!)
    
    % % %     im = im + 1;
    % % %     im2 = zeros(size(im) + [2 2 0]); % add borders
    % % %     im2(2:end-1, 2:end-1, :) = im; % put im at center
    im_black = (sum(im, 3) == 0);
    for nc = 1:3
        imc = im(:,:,nc);
        imc(im_black) = background_val(nc); % if black, then border value
        im(:,:,nc) = imc;
    end
    % % %     im = im2 - 1;
    
    fprintf('%g ', i)
    imwrite(im, [dir_out filesep file(1:end-4) num2str(i, '%03i') '.tif'], 'Compression', 'none')
end
fprintf('\n')

%%%

% for i = 1:N, im = imread('LN ing WT CD3.tif', i); im_black = (sum(im, 3) == 0); Nblack(i) = sum(im_black(:)); fprintf('%g ', i), end
