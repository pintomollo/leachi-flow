function new_names = equalize_stack(files, alpha)

%% function equalize_stack(file, alpha)
%
% Compute the mean intensity of each frame and equalise the stack
% so as to follow a parabolic variation I_fit(z) = a pi (r2 - z2)
% (2nd order Gauss Newton fit) as expected for a spheroid object.
%
% Correction modulated by the alpha exponent, for only partial correction:
% I_eq = I * (I_fit / mean(I))^alpha
%
% -- inputs --
% * file: image stack filename
% * alpha: coefficient for equalization (from 0 to 1, default 0.2, see equation above)
%
% -- output --
% Equalized images are automatically saved, with the extension '_filt.tif'.
%
% (c) Arnauld SERGE, 2011
%
% see also rendering_3D, filter_stack

  if nargin<1, files = '*.tif'; end

  new_names = {};
  dir_out = '_equalized';

  [files, out_path] = get_filenames(files, dir_out);

  Nk = length(files);
  if (Nk == 0), disp('nada??'), return, end
  %Nk = length(imfinfo(file));
  %kk = (1:Nk);

  new_names = files;

if nargin<2, alpha = 0.2; end

if (nargin <= 0)
  return;
end

im = imread(files{1}, 1);
nx = size(im, 1);
ny = size(im, 2);
nc = size(im, 3);
%nz = length(imfinfo(file));
nz = length(files);
mean_val = zeros(nz, nc);

fprintf(' Computing mean of image    ')
for iz = 1:nz
    im = imread(files{iz});
    for ic = 1:nc
        im_c = im(:, :, ic);
        mean_val(iz, ic) = mean(im_c(:));
    end
    fprintf('\b\b\b%3d',iz)
end
type = class(im);

order = 2;
if (mod(nz, 2)==0), frames = nz - 1; else frames = nz; end
I_fit = sgolayfilt(mean_val, order, frames);

% % fun = @(p,xdata) p(1)*xdata.^2+p(2)*xdata + p(3); lsqcurvefit(fun,x0,xdata,[0 0 0],ub)
I_fit2 = I_fit;
for ic = 1:nc
    if min(I_fit(:, ic))<0
        I_fit2(:, ic) = I_fit(:, ic) - min(I_fit(:, ic));
    end
end

subplot(121)
hold on
plot(fliplr(mean_val))
plot(fliplr(I_fit),':')
xlabel('z plane')
ylabel('mean fluorescent intensity')
title(['Equalizing ' file], 'interpreter', 'none')

im_eq = zeros(nx, ny, nc);
mean_eq = zeros(nz, nc);
%if isempty(dir('temp')), mkdir('temp'), end

fprintf('\n Equalizing image    ')
for iz = 1:nz
    filename = files{iz};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    im = double(imread(filename));
    for ic = 1:nc
        im_eq(:, :, ic) = im(:, :, ic) * real((I_fit2(iz, ic) / mean_val(iz, ic)).^alpha); % real for possible (negative values)^alpha generating irrationals..
        mean_eq(iz, ic) = mean(mean(im_eq(:, :, ic)));
    end
    
    %if iz == 1, imwrite(uint8(im_eq), [file(1:end-13) '_eq.tif'], 'tiff', 'Compression', 'none') % -13 to remove '_smoothed.tif'
    %else imwrite(uint8(im_eq), [file(1:end-13) '_eq.tif'], 'tiff', 'Compression','none', 'WriteMode', 'append')
    %end
    
    imwrite(cast(im_eq, type), new_name, 'TIFF');

    new_names{nk} = new_name;

    fprintf('\b\b\b%3d', iz)
end

plot(fliplr(mean_eq), '--')
a = axis; axis([a(1) a(2) 0 a(4)])
drawnow expose
pause(1)
%     saveas(gcf, ['3D' filesep file '_fluo_profil.fig'])

%%%
