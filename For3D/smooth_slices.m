function smooth_slices(file)

%% function smooth_slices(file)
%
% Compute & smooth object width & length over a zstack.
%
% -- input --
% file: image stack filename
%
% -- output --
% Smoothed images are automatically saved in a 'temp' subfolder, with
% the extension '_smoothed.tif'.
%
% see also rendering_3D
%
% (c) AS 2014


if isempty(dir('temp')), mkdir('temp'), end % dir_out = ['temp' filesep file(1:end-4) '_smooth']; if isempty(dir(dir_out)), mkdir(dir_out), end

if strfind(file, '_blue_brown.tif'), smooth_file = [file(1:end-15) '_smoothed.tif'];
else smooth_file = ['temp' filesep file(1:end-4) '_smoothed.tif'];
end

Nk = length(imfinfo(file));
kk = (1:Nk);

i_bary = zeros(Nk, 1);
j_bary = zeros(Nk, 1);
std_i = zeros(Nk, 1);
std_j = zeros(Nk, 1);
background_val = zeros(Nk, 3);

fprintf('filtering %i images for axial smoothing', Nk)

for nk = 1:Nk
    im = imread(file, nk);
    border_rows = squeeze([im(1,:,:) im(end,:,:)]);
    border_cols = squeeze([im(:,1,:); im(:,end,:)]);
    background_val(nk, :) = mean([border_rows; border_cols]);
    
    %% image pretreatment
    im3 = medfilt2(sum(im, 3)); % summing all channels to 'get all signal' (ponder color sum??)
    if mean(background_val(nk, :)) > mean(im(:)), im3 = max(im3(:)) - im3; end % invert if white bkground
    
    [Ni, Nj] = size(im3);
    ii = (1:Ni)';
    jj = (1:Nj);
    
    %% threshold
    threshold = opthr(im3(1:5:end), 0); % threshold = graythresh(im3)*max(im3(:));
    im3(im3 < threshold) = 0;
    % center and width computed for relevant values (above threshold <=> in the organ)
    val_i = mean(im3, 2);
    val_j = mean(im3);
    
    %% find typical size: barycenter & std (ponder by image values)
    i_bary(nk) = sum(val_i .* ii) / sum(val_i); % mean = ib = sum(Ii*i)/sum(Ii)
    j_bary(nk) = sum(val_j .* jj) / sum(val_j);
    std_i(nk) = sqrt(sum(val_i .* (ii - i_bary(nk)).^2) / sum(val_i)); % std = sqrt(sum(Ii * (i-ib)^2) / sum(Ii))
    std_j(nk) = sqrt(sum(val_j .* (jj - j_bary(nk)).^2) / sum(val_j));
    
    if mod(nk, 5) == 0 % show only some, to save time..
        imshow(uint8(im3/3)), title(['section ' num2str(nk) ', all channels summed']), hold on
        colormap(gray)
        plot(val_i, ii, 'r', jj, val_j, 'm', j_bary(nk), i_bary(nk), 'yo')
        
        xe = j_bary(nk)-std_j(nk);
        ye = i_bary(nk)-std_i(nk);
        we = 2*std_j(nk);
        he = 2*std_i(nk);
        rectangle('Position', [xe, ye, we, he], 'Curvature', [1, 1], 'edgecolor', 'y', 'linewidth', 2)% line([i_bary(nk)+std_i(nk) i_bary(nk)], [j_bary(nk) j_bary(nk)])
        legend({'mean intensity along i axis', 'mean intensity along j axis', 'elipse @ mean +/- std'})
        drawnow
    end
end

%% smooth shape descriptors => scaling factors
smooth_i_bary = smooth(i_bary); % using smooth, instead of expecting a spheroid shape
smooth_j_bary = smooth(j_bary);
smooth_std_i = smooth(std_i);
smooth_std_j = smooth(std_j);
i_scale = smooth_std_i ./ std_i; % scale = size_out / size_in
j_scale = smooth_std_j ./ std_j;

clf
subplot(2,3,1), plot([i_bary, smooth_i_bary]), title('i bary')
subplot(2,3,2), plot([j_bary, smooth_j_bary]), title('j bary')
subplot(2,3,4), plot([std_i, smooth_std_i]), title('std i')
subplot(2,3,5), plot([std_j, smooth_std_j]), title('std j')
subplot(2,3,3), plot(kk, i_scale*100, 'r', kk, j_scale*100, 'm', kk, ones(Nk)*100, 'k:')
legend({'corr i', 'corr j'})
figure(gcf), pause(.1)
warning off images:imshow:magnificationMustBeFitForDockedFigure

for nk = 1:Nk
    im = double(imread(file, nk));
    
    %% scale image by interp - image too large => pixels smaller & vice versa
    ii_corr = (ii - smooth_i_bary(nk)) / i_scale(nk) + i_bary(nk); % out : i_corr - <i_corr> = (i - <i>) / scale
    jj_corr = (jj - smooth_j_bary(nk)) / j_scale(nk) + j_bary(nk);
    [x_grid, y_grid] = meshgrid(jj_corr, ii_corr); % x = j & y = i!!
    im_corr = zeros(size(im)); % ([size(j_grid), 3]);
    
    for nc = 1:3
        im_corr_c = interp2(im(:, :, nc), x_grid, y_grid);
        im_corr_c(isnan(im_corr_c)) = background_val(nk, nc);
        im_corr(:, :, nc) = im_corr_c;
    end
    % subplot(232), imshow(im), subplot(235), imshow(im_corr)
    %     if mod(nk, 5) == 0
    subplot(2,3,6), imshow(uint8(im)), title(['section ' num2str(nk) ', resized & recentered'])
    hold on
    xe = j_bary(nk)-std_j(nk);
    ye = i_bary(nk)-std_i(nk);
    we = 2*std_j(nk);
    he = 2*std_i(nk);
    rectangle('Position', [xe, ye, we, he], 'Curvature', [1, 1], 'edgecolor', 'y', 'linewidth', 2)% line([i_bary(nk)+std_i(nk) i_bary(nk)], [j_bary(nk) j_bary(nk)])
    drawnow
    %     end
    
    %% save corrected image
    if nk == 1, imwrite(uint8(im_corr), smooth_file, 'tiff', 'Compression', 'none')
    else imwrite(uint8(im_corr), smooth_file, 'tiff', 'Compression', 'none', 'writemode', 'append')
    end %     file_out = [dir_out filesep file(1:end-4) '_smoothed' num2str(nk, '%03i') '.tif']; imwrite(uint8(im_corr), file_out, 'tiff', 'Compression', 'none')    %     bfsave(im,'myMultipageFile.tif', 'XYCZT')
end

%%%

% nk=1; Nk=1;
% im3=zeros(100); im3(40:60,40:60)=1;
% i_scale=2; j_scale=.5;
%
% im_corr3 = interp2(im3, j_grid, i_grid); imshow(im_corr3)