function filt_hsv(file)

%% function filt_hsv(file)
%
% Detect blue and brown stainings from HSV criteria
% and save them in blue & red channels of an "RGB" output.
%
% -- input --
% file: image stack filename
%
% -- output --
% Filtered & split images are automatically saved in a 'temp' subfolder,
% with the extension '_blue_brown.tif'.
%
% see also rendering_3D
%
% (c) AS 2014


Nk = length(imfinfo(file)); % or Nz...

warning off images:imshow:magnificationMustBeFitForDockedFigure

file_out = ['temp' filesep file(1:end-4) '_blue_brown.tif'];
if isempty(dir('temp')), mkdir('temp'), end

disp('Please close any Explorer Window that would display the current data folder, to avoid bug with writing images (indexing conflict..)')
fprintf('Splitting %i IHC images (blue, brown staining) in %s, according to HSV values:     ', Nk, file)

for nk = 1:Nk
    im = imread(file, nk);
    fprintf('\b\b\b%3i', nk)
    
    im = rgb2hsv(im);
    h = im(:, :, 1);
    s = im(:, :, 2);
    v = im(:, :, 3);
    
    %% brown
    s1 = s;
    v1 = v;
    
    ind = (v > 0.3); % & (s < 0.7)
    s1(ind) = 0;
    v1(ind) = 1;
    
    im2 = im;
    im2(:, :, 2) = s1;
    im2(:, :, 3) = v1;
    
    %% blue
    ind = (h < 0.5) | (h > 0.72) | (v < 0.3);
    s2 = s;
    v2 = v;
    s2(ind) = 0;
    v2(ind) = 1;
    
    im3 = im;
    im3(:, :, 2) = s2;
    im3(:, :, 3) = v2;
    
    %% output image
    im_out = zeros(size(im2));
    im_out(:,:,1) = 1-v1; % R
    im_out(:,:,2) = 1-v1; % G % duplicate for follicles
    im_out(:,:,3) = 1-v2; % B
    im_out = 256*im_out;
    
    %% plot
    %     if  (nk == 1) || (mod(nk, 10) == 0) % show only some?
    subplot(321), hist(h(1:100:end), 1000), title('H') % lines for thresholds...
    subplot(323), hist(s(1:100:end), 1000), title('S')
    subplot(325), hist(v(1:100:end), 1000), title('V')
    montage = [im; im2; im3];
    montage = [hsv2rgb(montage(1:10:end, 1:10:end, :)); im_out(1:10:end, 1:10:end, :)];
    subplot(122), imshow(flipdim(montage, 1))
    axis xy, ylabel('out / blue / brown / in'), title(['image  ' num2str(nk)])
    drawnow
    %     end
    
    %% save
    if (nk == 1), imwrite(uint8(im_out), file_out, 'tiff', 'Compression', 'none')
    else imwrite(uint8(im_out), file_out, 'tiff', 'Compression', 'none', 'writemode', 'append')
    end %     file_out = [dir_out filesep file(1:end-4) '_bb' num2str(nk, '%03i') '.tif'];    file_out = [dir_out filesep file(1:end-4) '_bb' num2str(nk, '%03i') '.tif'];    imwrite(uint8(256*im_out), file_out, 'tiff', 'Compression', 'none')
end % image
fprintf('\r')

%%%