function mystruct = get_struct(type, nstruct)
% GET_STRUCT retrieve custom data structures.
%
%   This function is designed as a centralization for the different 
%   complex structures used throughout the program so that they can
%   be edited consistently.
%
%   MYSTRUCT = GET_STRUCT(TYPE, SIZE) returns a matrix of size SIZE 
%   of the custom data structure of type TYPE. SIZE can be multi-
%   dimensional.
%
%   MYSTRUCT = GET_STRUCT(TYPE) returns one structure (SIZE = 1).
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 13.05.2014

  % Set the default size
  if (nargin == 1)
    nstruct = 1;
  end

  % Switch between all the different types
  switch type

    % The biological properties of the colonies (in um & s)
    case 'botrylloides_leachi'
      gmm = get_struct('gaussian_mixture');

      % The width distribution
      vessel = gmm;
      vessel.mu = [50.90; 88.27] / 2.5;
      vessel.sigma = cat(3, 84.06, 49.68) / 2.5;
      vessel.proportions = [0.934 0.066];

      % The radius / intensity of cells
      cells = gmm;
      cells.mu = [6.23 / 2.5 (255 - 175.01)];
      cells.sigma = [ 1.04 8.55; ...
                     8.55  1038.50] / 2.5;
      cells.proportions = 1;

      mystruct = struct('vessel_width', vessel, ...
                        'blood_cell', cells);

    % Structure used to parse the original files (reused in myrecording)
    case 'channel'
      mystruct = struct('color', 1, ...             % Color of the channel (index of 'colors')
                        'compression', 'none', ...  % Compression used for the temporary file
                        'cosmics', false, ...       % Remove the cosmic rays in the image (see imcosmics.m)
                        'detrend', false, ...       % Detrend the image (see imdetrend.m)
                        'file', '', ...             % Path to the original file
                        'fname', '', ...            % Name of the corresponding temporary file
                        'hot_pixels', true, ...     % Remove the hot pixels in the image (see imhotpixels.m)
                        'max', -Inf, ...            % Original maximum value used for rescaling
                        'min', Inf, ...             % Original minimum value used for rescaling
                        'metadata', '', ...         % Recordings metadata
                        'normalize', true, ...      % Normalize the whole stack
                        'type', 'dic');             % Type of channel

    % Structure storing the color codes used in the display
    case 'colors'
      mystruct = struct('colormaps', {{@gray, @redbluemap, @redgreencmap, @hot, @jet}}, ...     % The colors used for the images
                        'spots', {{'r','k','b','g','k'}}, ...                                   % The colors used for the detections
                        'spots_next', {{'b','r','k','b', 'w'}}, ...                             % The second colors for the detections
                        'status', {{'myg','mky','mkg','mbg','myg'}}, ...                        % The colors for the status of cells
                        'links', {{'y','k','k','y','y'}}, ...                                   % The colors for the links between cells
                        'paths', {{@redbluemap, @redgreencmap, @gray, @redbluemap, @gray}}, ... % The colors for the paths
                        'text', {{'r','k','b','g','k'}});                                       % The colors for the text in the movies

    % Structure to store detections from segmentations
    case 'detection'
      mystruct = struct('carth', NaN(1, 2), ...     % Cartesian position of the detections (Nx2)
                        'cluster', [], ...          % Temporal cluster containing the detection paths
                        'noise', [], ...            % Parameters of the image noise
                        'properties', []);          % Other properties computed on each detection, depending on the algorithm

    % Structure handling the export options
    case 'exporting'
      mystruct = struct('file_name', '', ...                % The name of the file to create
                        'low_duplicates', true, ...         % Do we use low duplicates paths ?
                        'full_cycles_only', false, ...      % Keep only tracks that both start and end with a division
                        'export_data', true, ...            % Do we export a CSV table of the data ?
                        'data_aligning_type', 'time', ...   % How do we align the paths ? (time/start/end)
                        'export_noise', false, ...          % Export the parameters of the noise in each image
                        'export_movie', false, ...          % Do we export an AVI movie ?
                        'movie_show_index', true, ...       % Do we display the track indexes in the movie ?
                        'movie_show_detection', true, ...   % Do we display the detected radii in the movie ?
                        'movie_show_paths', false, ...      % Do we display the connections of the tracking in the movie ?
                        'movie_show_reconstruction', false);% Do we display the reconstructed image using the detected spots ?

    case 'gaussian_mixture'
      mystruct = struct('mu', [], ...
                        'sigma', [], ...
                        'proportions', []);

    % The few parameters required to filter the image appropriately
    case 'image_filters'
      mystruct = struct('hot_pixels_threshold', 15, ...     % The threshold used to detect hot pixels using MEAN(pixels) +/- THRESH*STD(pixels)
                        'cosmic_rays_threshold', 15, ...    % The threshold used to detect cosmic rays as being separated from lower values by a gap larger than THRESH*MAD
                        'cosmic_rays_window_size', 10, ...  % The size of the overlapping tiles the algorithm works with
                        'detrend_meshpoints', 32);          % The number of positions over the image used for detrending

    % Structure containing the different parameters required for tracking spots
    case 'image_segmentation'
      mystruct = struct('filter_max_size', 10, ...         % Maximal radius (in um), see filter_spots.m
                        'filter_min_size', 2, ...          % Minimal radius (in um), see filter_spots.m
                        'filter_min_intensity', 15, ...    % Minimal intensity (multiplied by the variance of the noise), see filter_spots.m
                        'filter_overlap', 0.75, ...        % Minimal overlap (in percents) for spots to be fused together
                        'detrend_meshpoints', 32, ...      % The number of positions over the image used for detrending
                        'denoise_func', @gaussian_mex, ... % Function used to denoise the image (see imdenoise.m)
                        'denoise_size', -1,          ...   % Parameter used by the denoising function (-1: default value)
                        'denoise_remove_bkg', true, ...    % Removes the background uniform value (estimated using estimate_noise.m) ?
                        'atrous_max_size', 10, ...         % Maximal size of the spots to detect (in um), see imatrous.m
                        'atrous_thresh', 10, ...           % Threshold used to detect a valid spot as brighter than THRESH*MAD
                        'force_estimation', 1, ...         % Will force the estimation of the signal intensity on the raw data
                        'estimate_thresh', 1, ...          % Utilizes only the pixels brighter than this threshold (times noise variance) to estimate the spots
                        'estimate_niter', 15, ...          % Maximal number in the estimation procedure, see estimate_spots.m
                        'estimate_stop', 1e-2, ...         % Stopping criterion for the estimation procedure, see estimate_spots.m
                        'estimate_weight', 0.1, ...        % Convergence weight for the estiamtion procedure, see estimate_spots.m
                        'estimate_fit_position', false);   % Fit also the subpixel position of the spot (slower) ?

    case 'interpolation'
      mystruct = struct('segments', [], ...
                        'scaling', 1, ...
                        'speeds', [], ...
                        'dists', [], ...
                        'joints', [], ...
                        'parameters', []);

    case 'junction'
      mystruct = struct('threshold', 0, ...
                        'center', [], ...
                        'polygon', [], ...
                        'vector', []);

    case 'meshing'
      mystruct = struct('nodes', [], ...
                        'sorted', [], ...
                        'registraton', [], ...
                        'proportions', [], ...
                        'bounding_box', [], ...
                        'edges', []);

    % Structure used to handle the metadata provided by the microscope
    case 'metadata'
      mystruct = struct('acquisition_time', [], ... % Time of frame acquisition
                        'channels', {{}}, ...       % List of acquired channels
                        'exposure_time', [], ...    % Frame exposure time
                        'raw_data', '', ...         % Raw metadata string
                        'z_position', []);          % Z-position of the frame

    case 'MicroMos'
      mystruct = struct('ImageFolder', '', ... % absolute or relative path where the images are stored inside the computer.
                        'ImageBaseName', '', ... % name of the images to be stitched, without the final cardinal number.
                        'ImageIndexs', [], ... % arabic cardinal numbers of the specific images to be stitched between the ones present in the "ImageFolder". The numbers must be reported in the order of the images to be stitched.
                        'flag_Color', 1, ... % (by default: 1). 0 to obtain the final mosaic in grey levels (also if the original images are RGB). 1 to obtain the final mosaic in RGB (only if the oiginal images are RGB). 
                        'flag_SeekBestImages', 0, ... % (by default: 1). The input frame are pre-selected to optimize the registration task. If this parameter is set to 1, not all the images pointed out are stitched in the final mosaic.
                        'flag_BleachingCorrection', 0, ... % (by default: 0 for bright field and phase contrast images, 1 for fluorescent images). To correct intensity decay, due to the photo-bleaching effect, between the images to be stitched. This phenomenon happen especially if the images to be stitched are fluorescent images. 0 = no bleaching correction; 1 = yes bleaching correction.
                        'ShiftEstimationMode', 2, ... % 1: PhaseCorrelation, 2: CornerClustering, 0: none
                        'RegistrationMode', 2, ... % (Piccinini, 2013: proj and aff identical, trans better if very little overlap); (by default: uint8(2)). To choose the registration model that must be used to register the images. The registration model can be chosen between projective (suggested), affine and translative. Set: uint8(0) = 'projective' or uint8(1) = 'affine' or uint8(2) = 'translative' (computed at sub-pixel accuracy using the Lukas-Kanade feature tracker) or uint8(3) = 'translative' (computed at pixel accuracy using the phase-correlation algorithm only).
                        'flag_Blending', 1, ... % (Piccinini, 2013: Not required for widefield, not exactly correct, looks much better with biquadratic...) 0 = no blending; 1 = blending using a biquadratic function with the highest value in the image's centre (seams are only attenuated). 2 = linear blending with boundary values to completely avoid seams in the stitching regions (slow computation using pixel-based interpolation). 3 = linear blending with boundary values to completely avoid seams in the stitching regions (fast computation using Delaunay triangulation).
                        'flag_GriddedAcquisition', 0, ... % 1 if the mosaic is a predictible NxM grid
                        'GridMode', 0, ... % 0: Row by row, from left to right
                        'flag_WhiteBalancing', 1, ... % (??) 1 to perform the white balancing of the output mosaic using the mosaic itself as colour reference. 2 to perform the white balancing of the output mosaic loading an external 3-channel image (a RGB image) that must be copied in the folder called "WHITEBALANCING". 0 otherwise.
                        'flag_FlatField', 1, ... % (Piccinini, 2013: crucial) 1 to flat-field correct the images using an input vignetting function. The vignetting function must be saved as matlab matrix in the folder named: "VIGNETTINGFUNCTION". In the "VIGNETTINGFUNCTION" folder must contain at maximum one vignetting function file.
                        'flag_FrameToMosaic', 1, ... % (Piccinini, 2013: F2M clearly better) (by default: 1). 0 for registering the images according to the Frame-to-Frame registration approach; 1 (suggested) for registering the images according to the Frame-to-Mosaic registration approach.
                        'RANSACerror', 2, ... % maximum subpixel reprojection error used inside RANSAC algorithm. We suggest 2 pixels.
                        'flag_PhaseCorrelationOnly', 0, ... % It can assume values 0 or 1. 1 means that the images are registered according to the Phase Correlation ALgorithm only.
                        'flag_PCglobalORlocal', 0, ... % It can assume values 0 or 1. 0 means that the metric used to determine the best shift inside the Phase Correlation ALgorithm is the global RMSE performed on the whole overlapping region. 1 means that the used metric is the RMSE performed using only the pixels with highest value. 
                        'PCscaleFactor', 1, ... % to speed up the computational processes. Image rescale factor applyed to optionally resize the images. The value must be a positive integer. E.g.: ceil(2) to obtain half-sized images than the original images. 1 (suggested) means: rescaling not active.
                        'flag_ComputeRegistrations', 1, ... % 0 for loading an external registration matrix to stitch the images. In case, the registration matrix must be saved as 3x3xn (n=number of images to be registered) Matlab matrix in the folder named: "REGISTRATIONMATRIX".
                        'InterpolationMode', 'bilinear', ... % interpolation used to warp the image to be stitched. It can be: 'bicubic' or 'bilinear' (suggested) or 'nearest'.
                        'flag_LookUpTable', 1, ... % to use a 256 levels Look-Up-Table to map the grey levels into a single RGB conversion.
                        'ScaleFactor', 1, ... % to speed up the computational processes. Image rescale factor applyed to optionally resize the images. The value must be a positive integer. E.g.: ceil(2) to obtain half-sized images than the original images. 1 (suggested) means: rescaling not active.
                        'PixelAccuracy', 0); % (by default: 0). Pixel-accuracy of the output mosaic. 0 = double precision (slower computation, higher precision), 1 = single precision (faster computation, lower precision). If your computer have not enough memory to build the mosaic, we suggest to try to set this parameter at the value 1 (single precision).

    % Global structure of a recording and analysis
    case 'myrecording'
      mychannel = get_struct('channel', 0);
      mysegment = get_struct('segmentation', 0);
      mytracks = get_struct('tracking', 0);
      mystruct = struct('channels', mychannel, ...      % Channels of the recording
                        'segmentations', mysegment, ... % Segmentation data
                        'trackings', mytracks, ...      % Tracking data
                        'experiment', '');              % Name of the experiment

    % Global structure of options/parameters for an analysis
    case 'options'
      myfilt = get_struct('image_filters');
      mysegm = get_struct('image_segmentation');
      mytrac = get_struct('spot_tracking');
      mytrkf = get_struct('tracks_filtering');
      mystruct = struct('config_files', {{}}, ...       % The various configuration files loaded
                        'binning', 1, ...               % Pixel binning used during acquisition
                        'ccd_pixel_size', 16, ...       % X-Y size of the pixels in um (of the CCD camera, without magnification)
                        'magnification', 20, ...        % Magnification of the objective of the microscope
                        'spot_tracking', mytrac, ...    % Parameters for tracking the spots
                        'filtering', myfilt, ...        % Parameters for filtering the recordings
                        'tracks_filtering', mytrkf, ... % Parameters for filtering the tracks
                        'pixel_size', -1, ...           % X-Y size of the pixels in um (computed as ccd_pixel_size / magnification)
                        'segmenting', mysegm, ...       % Parameters for segmenting the recordings
                        'time_interval', 300, ...       % Time interval between frames (in seconds)
                        'verbosity', 2);                % Verbosity level of the analysis

    % Structure used to segment a channel
    case 'segmentation'
      mydetection = get_struct('detection',0);
      mystruct = struct('denoise', true, ...           % Denoise the segmentation (see imdenoise) ?
                        'detrend', false, ...          % Detrend the segmentation (see imdetrend.m) ?
                        'filter_spots', true, ...      % Filter the spots (see filter_spots.m) ?
                        'detections', mydetection, ... % The structure used to store the resulting detections
                        'type', {{}});                 % The type of segmentation

    case 'simulation'
      mystruct = struct('image_size', [512 512], ...
                        'image_noise', 0.05, ...
                        'dt', 1, ...
                        'dmax', 0.5, ...
                        'init_simulation', @create_vessel, ...
                        'init_params', [1 2], ...
                        'create_cells', @vessel_init, ...
                        'creation_params', [], ...
                        'move_cells', @vascular_movement, ...
                        'movement_params', 200, ...
                        'remesh_cells', @threshold_remesh, ...
                        'remesh_params', 0.75, ...
                        'ndims', 2, ...
                        'nprops', 2, ...
                        'pixel_size', 5, ....
                        'duration', 300, ...
                        'outside_ridge', 0.4, ...
                        'cell_speed', 20, ...
                        'cell_size', 2, ...
                        'cell_variation', 0.1, ...
                        'cell_density', 0.0005);

    % Structure containing the different parameters required for tracking spots
    case 'spot_tracking'
      mystruct = struct('spot_max_speed', 0.05, ...           % Maximal speed of displacement of a spot (in um/s)
                        'allow_branching_gap', false, ...     % Allows merging/splitting to occur at the same time as gap closing ?
                        'min_section_length', 5, ...          % The minimum number of contiguous frames a path "section" should last to be kept for merging/spliting/gap closing
                        'bridging_max_gap', 3, ...            % Considered number of frames for the gap closing algorithm (see track_spots.m)
                        'max_intensity_ratio', Inf, ...       % Defines an upper bound to the allowed signal ratios (see track_spots.m)
                        'bridging_max_dist', Inf, ...         % Maximal distance throughout the gap (see track_spots.m)
                        'bridging_function', @bridging_cost_sparse_mex, ...   % Function used to measure the gap-closing weight
                        'joining_function', @joining_cost_sparse_mex, ...     % Function used to measure the joinging weight
                        'splitting_function', @splitting_cost_sparse_mex, ... % Function used to measure the splitting weight
                        'linking_function', @linking_cost_sparse_mex); ...    % Function used to measure the frame-to-frame linking

    % The trackings as stored after segmentation
    case 'tracking'
      mydetection = get_struct('detection',0);
      mystruct = struct('reestimate_spots', true, ...      % Do we reestimate the newly interpolated spots ?
                        'filtered', mydetection, ...       % The structure used to store the detections after filtering
                        'detections', mydetection);        % The structure used to store the resulting detections

                        % Functions that are not implemented yet :
                        %'force_cell_behavior', true, ...   % Prevent fusion and appearance of spots
                        %'post_processing_funcs', {{}}, ... % Allow to post-process paths

    % The options for filtering tracks
    case 'tracks_filtering'
      mystruct = struct('interpolate', true, ...        % Interpolate the position of cells to fill gaps ?
                        'max_zip_length', 3, ...        % A "zip" is defined as a cell splitting and merging back together. This defines the maximum number of frames this separation can last to be closed
                        'min_path_length', 10);         % The minimum number of frames a path should last to be kept

    case 'vessel'
      mystruct = struct('center', [], ...
                        'border', [], ...
                        'junction', [], ...
                        'property', [], ...
                        'mesh', []);

    % If the required type of structure has not been implemented, return an empty one
    otherwise
      mystruct = struct();
  end

  % Compute the pixel size
  mystruct = set_pixel_size(mystruct);

  % Repeat the structure to fit the size (nstruct can be multi-dimensional)
  mystruct = repmat(mystruct, nstruct);

  return;
end
