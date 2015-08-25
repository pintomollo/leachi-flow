function rendering_3D(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function rendering_3D(params)
%
%% This function is the main one for Full Organ Rendering in 3D (For3D)
%% from a complete collection of immunohistologically labelled slices
%
%% Preliminary steps, in ImageJ
% 1/ adjust contrast & remove background for each color if necessary
% 2/ correct/remove artifacts if necessary (flipped sections, fat tissue, etc.)
% 3/ run the StackReg plugin (Thevenaz et al, EPFL)
% 4/ crop & save as multi-tif
%
% Color thresholds may be assessed in ImageJ from montage image
% (=> seing all images of stack at once)
%
%% -- input --
% params: parameters can be set either using the GUI or directly by the command line:
% i.e.
% rendering_3D({'*.tif'; '0'; '[130 -1 -1]'})
% => testing default threshold values for green and blue channels
%
% rendering_3D({'*.tif'; '0'; '[130 35 200]'})
% => using manually set thresholds
%
% default:
% rendering_3D({'*.tif';'0';'[-1 -1 -1]';'6.5';'22';'0.2';'[1/3 1/20 1/20]';'[0 0]';'0'})
%
%% -- output --
% For an input stack named 'input_file.tif', 3D images (as well as intermediate steps)
% are automatically saved in a subfolder named'temp/input_file_3D', with the extension
% '_3D_NN.png' (NN ranging from 0 to 45).
%
% see also filter_stack, equalize_stack, FOR3D_dialog_box
%
%% (c) Arnauld SERGE, 2011
% Contact: arnauld.serge@univ-mrs.fr
%
%% See ref:
%% J Immunol. 2013 Jan 15;190(2):586-96. doi: 10.4049/jimmunol.1200119. Epub 2012 Dec 17
%
%% Three-dimensional visualization of the mouse thymus organization in health and immunodeficiency.
%% Irla M, Guenot J, Sealy G, Reith W, Imhof BA, Sergé A.
%% Department of Pathology and Immunology, Faculty of Medicine, University of Geneva, Switzerland
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
  params = get_struct('For3D');
end

if isempty(params.filename), return, end % cancel by user

%Nsampling = floor(params.slice_width/params.pixel_size/1.5); % assuming ideally dz = 1.5dx for resoluting voxel after sampling.
%if Nsampling < 1, Nsampling = 1; end % Thymus (screen capture at 2x): floor(22/6.5/1.5) = 2, Thymus (3x): floor(22/4.33/1.5) = 3, LN (8x): floor(22/2.6/1.5) = 5

%{
%% select files
files = dir(filename);
Nf = length(dir(filename));

do_filt = isempty(strfind(filename, '_filt.tif')); % file has to be, or already filtered?

if do_filt
    ok = ones(1, Nf);
    for nf = 1:Nf % sort files: keep only filt files
        if strfind(files(nf).name, '_filt.tif'), ok(nf) = 0; end
    end
    files = files(ok==1);
    Nf = sum(ok);
end

if Nf == 0, disp(['Sorry, ' filename ' not found in ' cd]), end
thresholds = repmat(thresholds, Nf, 1); % allowing giving different values / file?
%}

  files = params.filename;
  %new_names = {};
  %dir_out = '_resampled';

  %[files, out_path] = get_filenames(files, dir_out);
  [files] = get_filenames(files);

  Nf = length(files);
  if (Nf == 0), disp('nada??'), return, end

  thresholds = params.thresholds;
  %new_names = files;

  im = imread(files{1});
  [Nx, Ny, Nz] = size(im);
  Nc = Nf;
  nf = 1;

%% loop/files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for nf = 1:Nf
    
    %{
    file = files(nf).name;
    if do_filt, file_out = [file(1:end-4) '_3D'];
    else file_out = [file(1:end-9) '_3D'];
    end
    dir_out = ['temp' filesep file_out];
    
    if ~isempty(dir([dir_out filesep file_out '45.png']))
        disp([file ' already done'])
        continue
    end
    
    fig_file = ['temp' filesep file_out '.fig'];
    % % % % % %     if isempty(dir(fig_file)) % work with already saved fig?
    
    if nf == 1
        fprintf('pixel size: %.2f um before and %.2f um after %i times sampling\n', pixel_size, pixel_size*Nsampling, Nsampling)
    end
    
    if do_filt
        if strfind(file, '_filt.tif'), continue, end % filtered files..
        
        file_blue_brown = ['temp' filesep file(1:end-4) '_blue_brown.tif'];
        smooth_file = ['temp' filesep file(1:end-4) '_smoothed.tif'];
        eq_file = ['temp' filesep file(1:end-4) '_eq.tif'];
        filt_file = [file(1:end-4) '_filt.tif'];
        
        if isempty(dir(smooth_file)) || isempty(dir(eq_file)) || isempty(dir(filt_file))
            figure('windowstyle', 'docked');
        end
        
        %% filter hsv staining
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if detect_IHC
            filt_hsv(file) % generate blue_brown stack
            unsmoothed_file = file_blue_brown;
        else
            unsmoothed_file = file;
        end
        
        %% smooth section sizes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(dir(smooth_file))
            smooth_slices(unsmoothed_file) % generate smoothed stack
        end
        
        %% equalize & filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stk = filter_stack(file, alpha, Nsampling, border_erosion/pixel_size); % generate eq & filt stack
        if isempty(stk), return, end % incompatiblity with old filt file?
        
    else % starting directly with filt file
        im = imread(file, 'Index', 1);
        [Nx, Ny, Nc] = size(im);
        Nz = length(imfinfo(file));
        stk = zeros(Nx, Ny, Nz, Nc);
        for nz = 1:Nz, stk(:, :, nz, :) = imread(file, 'Index', nz); end
    end
    
    [Nx_out, Ny_out, Nz_out, Nc] = size(stk);
    %}
    
    % % % % %     if detect_IHC
    % % % % %         %% extra median filter for B zones
    % % % % %         Msize = [7 7]; % odd size is better
    % % % % %         fprintf('\n Extra median filtering, over %ix%i pixels, %s, image    ', Msize, file)
    % % % % %         for nz = 1:Nz_out
    % % % % %             stk(:, :, nz, 2) = medfilt2(stk(:, :, nz, 2), Msize);
    % % % % %             fprintf('\b\b\b%3d', nz)
    % % % % %         end
    % % % % %         fprintf('\n')
    % % % % %     end
    
    volume_c = zeros(nf, nc);

    %% thresholds auto
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nc = 1:Nc
        [stk, junk, type] = load_sparse_stack(files{nc}, params.sparse_thresholds(nc));

        if thresholds(nf, nc) == -1
            fprintf('finding automatic threshold %g. Iteration...', nc)

            T1 = mygraythresh(stk, type);
            T2 = full(opthr(reshape(stk, Nx, [])));
            %vals = stk(:, :, :, nc);
            %T1 = graythresh(vals)*max(vals(:));
            %T2 = opthr(reshape(vals, Nx_out, Ny_out*Nz_out));
            
            if Nc>1
                switch nc % tricky combination, might be optimized for a given dataset, or set manually
                    case 1 % red
                        if detect_IHC
                            a1 = 0.3; a2 =0.7;
                        else % intermediate threshold to get only medulla, avoiding non spe cortex
                            a1 = 0.7; a2 =0.3;
                        end
                        T = min((a1*T1 + a2*T2), 255);
                        fprintf('Red threshold computation: T1 = %.1f, T2 = %.1f, %g*T1 + %g*T2 = %.1f\r', T1, T2, a1, a2, T)
                        thresholds(nf, nc) = T;
                        
                    case 2 % green 
                        if detect_IHC
                            a1 = 0.3; a2 =0.7; a3 = 2;
                            T = min((a1*T1 + a2*T2)*a3, 255);
                            fprintf('Green threshold computation:  T1 = %.1f, T2 = %.1f, (%g*T1 + %g*T2)*%g = %.1f\r', T1, T2, a1, a2, a3, T)
                        else % low threshold to get all staining
                            T = T2; %/2;
                            fprintf('Green threshold computation:  T = %.1f\r', T)
                        end
                        thresholds(nf, nc) = T;
                        
                    case 3 % blue 
                        if detect_IHC
                            T = T2;
                            fprintf('Blue threshold computation:  T = %.1f\r', T2)
                        else % high threshold to get only capsule
                            a1 = 0.8; a2 =0.2; a3 = 1.5;
                            T = min((a1*T1 + a2*T2)*a3, 255); % (0.8*T1 + 0.2*T2)*1.5;
                            fprintf('Blue threshold computation:  T1 = %.1f, T2 = %.1f, (%g*T1 + %g*T2)*%g = %.1f\r', T1, T2, a1, a2, a3, T)
                        end
                        thresholds(nf, nc) = T;
                end
            else % gray, Nc = 1
                thresholds(nf) = mean([T1 T2]);
            end
            
            %% plot thresholds
            col = 'rgb';
            vol = zeros(255, 1);
            stk = stk(:);
            for th = 1:255, vol(th) = sum(stk<th); end
            volume_c(nc) = sum(stk > thresholds(nf, nc));
            
            th = round(thresholds(nf, nc));
            if th>0
                subplot(222)
                semilogy(vol, 'color', col(nc))
                hold on, axis tight
                semilogy([th th], [min(vol) vol(th)], 'color', col(nc))
                xlabel('threshold value'), ylabel('volume')
                title(['Testing threshold for ' file], 'interpreter', 'none')
                
                subplot(224)
                diffv = diff(vol);
                semilogy(diffv, 'color', col(nc))
                hold on, axis tight
                semilogy([th th], [1 diffv(th)], 'color', col(nc))
                xlabel('threshold value'), ylabel('derivate')
                drawnow expose
                pause(1)
            end
        end % if thresholds(nf, nc) == -1
    end % for nc = 1:Nc
    
    %voxel_size = (pixel_size*Nsampling)^2*slice_width*1e-9; % convert from um3 to mm3
    voxel_size = (pixel_size)^2*slice_width*1e-9; % convert from um3 to mm3
    for nc = 1:Nc
        %vals = stk(:, :, :, nc);
        %volume_c = sum(vals(:) > thresholds(nf, nc));
        fprintf('Channel %i: threshold = %.0f, volume = %i voxels = %.2f mm3\n', nc, thresholds(nf, nc), volume_c(nc), volume_c(nc)*voxel_size)
    end
    
    % % % % %     if detect_IHC % extra smoothing?
    % % % % %         [xx, yy, zz] = ndgrid(-5:5);
    % % % % %         nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 5; sphere r = 5, v = 524
    % % % % %         stk(:, :, :, 1) = imdilate(imerode(stk(:, :, :, 1), nhood), nhood);
    % % % % %     end
    
    %% ****** generate 3D view ******
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('computing 3D...\n')
    [xx, yy, zz] = meshgrid(1:Ny_out, 1:Nx_out, 1:Nz_out);
    
    h = figure('Color', 'w', 'Position', [80 80 560*2 420*2]); % for larger frames
    
    p = zeros(1, Nc);
    clr = 'rgb';
    
    %%  plot isosurf 
    for nc = 1:Nc
        fprintf('generating surface for channel %i...\n', nc)

        [stk, junk, type] = load_sparse_stack(files{nc}, params.sparse_thresholds(nc));

        %vals = squeeze(stk(:, :, :, nc));
        p(nc) = patch(isosurface(xx, yy, zz, stk, thresholds(nf, nc), 'verbose'));
        fprintf('rendering surface for channel %i...\n', nc)
        isonormals(xx, yy, zz, stk, p(nc))
        if size(clr, 1)>1
            set(p(nc), 'FaceColor', clr(nc, :), 'EdgeColor', 'none', 'FaceAlpha', transparency(nc))
        else
            set(p(nc), 'FaceColor', clr(nc), 'EdgeColor', 'none', 'FaceAlpha', transparency(nc))
        end
    end
    clear xx yy zz stk
    
    %% 3D volume
    view(3); axis tight, grid on
    camlight, lighting gouraud
    
    %% axes
    %xy_calib = pixel_size*Nsampling; % um
    xy_calib = pixel_size; % um
    z_calib = slice_width;
    daspect([1/xy_calib 1/xy_calib 1/z_calib])
    if (pixel_size > 1), axe_unit = 1000; axe_calib = 1; % um => 1 mm
    elseif detect_IHC, axe_unit = 100; axe_calib = 10; % 100 um, grad every 10 lines, hence every mm
    elseif strfind(filename, 'testis'), axe_unit = 1000; axe_calib = 10; % 1 mm, grad every 10 lines, hence every 10 mm
    else axe_unit = 1; axe_calib = 1;
    end
    
    a = axis;
    xmax = ceil(a(2)*xy_calib/axe_unit);
    ymax = ceil(a(4)*xy_calib/axe_unit);
    zmax = ceil(a(6)*z_calib/axe_unit);
    axis([0 xmax/xy_calib 0 ymax/xy_calib 0 zmax/z_calib]*axe_unit)
    % xmax = floor(a(2)*xy_calib/axe_unit); ymax = floor(a(4)*xy_calib/axe_unit); zmax = floor(a(6)*z_calib/axe_unit);
    set(gca, 'XTick', (0:xmax)*axe_unit/xy_calib, 'YTick', (0:ymax)*axe_unit/xy_calib, 'ZTick', (0:zmax)*axe_unit/z_calib)
    
    Xticklbl = cell(xmax+1, 1); % initially empty
    Yticklbl = cell(ymax+1, 1);
    Zticklbl = cell(zmax+1, 1);
    Xticklbl(1:axe_calib:end) = num2cell((0:axe_calib:xmax)/axe_calib); % use only one every axe_calib
    Yticklbl(1:axe_calib:end) = num2cell((0:axe_calib:ymax)/axe_calib);
    Zticklbl(1:axe_calib:end) = num2cell((0:axe_calib:zmax)/axe_calib);
    
    set(gca, 'XTicklabel', Xticklbl, 'YTicklabel', Yticklbl, 'ZTicklabel', Zticklbl)
    set(gca, 'FontSize', 16)
    
    if strfind(file, '_no_corner'), plot_3D_corner(Nx, Ny, Nz), end
    
    saveas(h, fig_file) % axis off, saveas(h, [file_out '.png']), axis on
    
    % % % % % %     else
    % % % % % %         h = load(fig_file);
    
    %% ****** video ******
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % using saveas (or imwrite, writeVideo, 40% faster but requires no screensaver, no other window on top..)
    az = 0; el = 90;
    if isempty(dir(dir_out)), mkdir(dir_out), end
    fprintf('saving view   '), tic % Video3D = VideoWriter([file '.avi']); Video3D.FrameRate = 7; open(Video3D);
    for i = 0:45
        view(az-2*i, el-2*i)
        % %         frame = getframe(h); writeVideo(Video3D, frame);
        % %         if i == 0, imwrite(uint8(frame.cdata), file_out, 'tiff', 'Compression', 'none')
        % %         else imwrite(uint8(frame.cdata), file_out, 'tiff', 'Compression', 'none', 'writemode', 'append'), end
        saveas(h, [dir_out filesep file_out num2str(i, '%02i') '.png']);
        fprintf('\b\b\b%3d', i)
    end
    fprintf('\n'), toc % close(Video3D);
    % % % % % %     end %     if isempty(dir(fig_file))
%end % /files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%

% adjusting colors a posteriori (ad libido!) i.e. with these commands:
%
% az = 0; el = 90; i = 15; view(az-2*i, el-2*i)
% h = get(gca, 'Children'); % => 'light'    'patch'    'patch'
% with no corner => 10 vals: light, 6 lines, 3 patches
% set(h(2), 'facecolor', [1 0 .5]), set(h(3), 'facealpha', 0.3)
% set(h(2), 'facecolor', [1 .5 0]), set(h(2), 'facealpha', 0.2), set(h(3), 'facecolor', [0 0 1]), set(h(3), 'facealpha', 0.5)
% saveas(gcf, 'pretty colors in 3D.png')

%%%
