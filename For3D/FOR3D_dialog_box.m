function [filename, detect_IHC, thresholds, pixel_size, slice_width,...
    alpha, transparency, border_erosion] = FOR3D_dialog_box(params)

%% function [filename, detect_IHC, thresholds, pixel_size, slice_width,...
%%     alpha, transparency, folder_3D, border_erosion] = FOR3D_dialog_box(params)
%
% User interface for starting rendering_3D
%
% -- inputs --
% params: list of For3D parameters (default: set by user, through the GUI)(0: use default values)
%
% -- outputs --
% [filename ... border_erosion]: list of For3D parameters, eventually modified through the GUI
% default: params = {'*.tif';'0';'[-1 -1 -1]';'6.5';'22';'0.2';'[1/3 1/20 1/20]';'[0 0]';'0'};
%
% (c) Arnauld SERGE, 2011
%
% see also rendering_3D


name_parameters = {...
    'Filename (default: all tif files)';...
    'Detect IHC staining (1=yes 0=no)';...
    'Red Green Blue thresholds (-1=auto)';...
    'Pixel size (um)';...
    'Section thickness (um)';...
    'Equalization factor (0=none, 1=full)';...
    'Red Green Blue transparency ratio';... % use [1/3 1/3 1/10] for LN, for same transparency in red and green
    'Red Green border erosion (um, ref=Blue)';...
    };


%% **** IF RECQUIRED, EDIT THESE LINES TO CHANGE DEFAULT VALUES ****
params_default = {...
    '*.tif';...
	'0';...
    '[-1 -1 -1]';...
    '6.5';...
    '22';...
    '0.2';...
    '[1/3 1/20 1/20]';...
    '[0 0]';...
    };

%% *******************

if isempty(params)
    params = inputdlg(name_parameters, 'For3D Parameters', 1, params_default);
    if isempty(params)
        filename = ''; detect_IHC = ''; thresholds = ''; pixel_size = ''; slice_width = ''; alpha = ''; transparency = ''; border_erosion = '';
        return
    end
else
    if isnumeric(params) % params = {...} format cell par defaut, numerique (1,2...) pour raccourcis...
        if (params == 0)
            params = params_default;
        elseif (params == 1)
            params = {'EG7 2RB no eros_no_corner_filt.tif';'0';'[12 0 12]';'6.5';'22';'0.5';'[1/3 0 1/16]';'[0 0]';'1'}; % tumor
        elseif (params == 2)
%             params = {'*.tif';'0';'[-1 -1 -1]';'2.55';'33';'0.2';'[1/3 1/20 1/20]';'[11 0]'}; % thymus
            params = {'*no_corner_filt.tif';'0';'[82 40 160]';'4.33';'22';'0.2';'[1/3 1/20 1/20]';'[11 0]'};
        elseif (params == 3)
            params = {'*.tif';'1';'[200 25 60]';'0.9';'5.5';'0.2';'[.2 .1 .5]';'[0 0]'}; % LN B220
            % rendering_3D({'*CD3.tif';'1';'[45 70 -1]';'0.9';'5.5';'0.2';'[.2 .1 .5]';'[0 0]'}) % LN CD3
        end
    end
    if length(params)<length(params_default) % if only first values given, complete with default
        params = [params(:); params_default(length(params)+1:end)];
    end
end

filename = params{1};
detect_IHC = eval(params{2});
thresholds = eval(params{3});% thresholds(strcmpi('auto', thresholds)) = {'-1'};% thresholds(strcmpi('saved', thresholds)) = {'-2'};% thresholds = str2double(thresholds);
pixel_size = eval(params{4});
slice_width = eval(params{5});
alpha = eval(params{6});
transparency = eval(params{7}); % array of 3 values
% % redo = str2double(params{9});% % do_plot = str2double(params{10});% % save_3D = str2double(params{11}); folder_3D = params{8}; % string!
border_erosion = eval(params{8});

if isempty(filename), disp('cancelled'), return, end % cancel by user

save_params(name_parameters, params)

%%%