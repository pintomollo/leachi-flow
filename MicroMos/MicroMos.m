function [FileName] = MicroMos(varargin)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 July 2013
% NAME: MicroMos
% 
% To build the mosaic according to the input parameters.
%
% INPUT:
%  parameters       Structure containing the parameters set into the files
%                   "ParametersByUser.m" and "ParametersDefault.m". Please 
%                   read carefully the two files to understand which 
%                   parameters are required.
%
% OUTPUT:
%  Mosaic           Matrix of the mosaic built stitching the "n" input 
%                   images according to the input parameters.
%  MaskOverlap      Matrix of the same size of "Mosaic". Each pixel reports
%                   the index of the last image (between the "n" images
%                   registered into the mosaic) that wrote in that specific
%                   position.
%  MatricesGLOBAL   3x3xn-1 ("n" = number of images composing the mosaic) 
%                   registration matrices used to warp and stitch into the 
%                   mosaic the last n-1 images.
%
% See also CheckPointsAndNAN, ColinearDegeneration, FlatFieldCorrection,
% ImagesForFTM, Metric_MSE, Metric_RMSE, Metric_SNR, Metric_UQI, MicroMos,
% MosaicAssessment, MosaicUpdating, ParametersByUser, ParametersDefault,
% PointsDetectionTracking, RansacAffine, RANSACmodified, RansacProjective,
% RansacTranslative, RegistrationErrorAssessment, ReprojectionError,
% ShiftByPhaseCorrelation, SlideShow, START, WarpingRegistrationMode

% CVG (Computer Vision Group) Toolbox
% Copyright © 2012 Filippo Piccinini, Alessandro Bevilacqua, 
% Advanced Research Center on Electronic Systems (ARCES), 
% University of Bologna, Italy. All rights reserved.
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License version 2 (or higher) 
% as published by the Free Software Foundation. This program is 
% distributed WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
% General Public License for more details.
%
% Massive speed, memory and ease of use improvements
% Simon Blanchoud, Wilson lab, University of Otago, 2015

  % Inputs processing
  if (length(varargin) > 0 && isstruct(varargin{1}))
    parameters = update_structure(varargin{1}, 'MicroMos');
    varargin(1) = [];
  else
    parameters = get_struct('MicroMos');
  end

  % Now we check that the parameters were provided in pairs
  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error 'Properties pairs must come in PAIRS.';
  end

  % Loop over the pairs of parameters
  for i = 1:npairs
    % If the parameter exists in opts we simply assign it the
    % provided value
    if (isfield(parameters, varargin{2*i - 1}))
      parameters.(varargin{2*i - 1}) = varargin{2*i};

    % Or ignore it !
    else
      warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
    end
  end

  %% PARAMETERS SETTTING
  GLOBAL = eye(3,3);
  MatricesGLOBAL = GLOBAL;
  LookUpTable = NaN.*ones(4,256);

  if parameters.PixelAccuracy == 0
      fPixelAccuracy = @double;
  else
      fPixelAccuracy = @single;
  end

  if (isempty(parameters.ImageFolder))
    parameters.ImageFolder = uigetdir('Movies','Select the directory containing the mosaic');
  end

  FileName = {};
  if (isempty(parameters.ImageBaseName))
    files = dir(fullfile(parameters.ImageFolder, '*'));

    fname = '';
    for i=1:length(files)
      if (files(i).name(1)~='.')
        if (files(i).isdir)

          others = dir(fullfile(parameters.ImageFolder, [files(i).name '.*']));

          if (isempty(others))
            tmp_params = parameters;
            tmp_params.ImageFolder = fullfile(tmp_params.ImageFolder, files(i).name);
            disp(['Performing mosaicing on subfolder ''' files(i).name '''']);

            try
              new_file = MicroMos(tmp_params);
              if (~isempty(new_file))
                FileName{end+1} = new_file;
              end
            catch
              err = lasterror();
              warning(['Error during the analysis:']);
              print_all(err);
            end
          end
        else
          [file_path, file_name, file_ext] = fileparts(files(i).name);
          if (exist(fullfile(parameters.ImageFolder, file_name), 'dir'))
            FileName{end+1} = fullfile(parameters.ImageFolder, files(i).name);
          else
            if (isempty(fname))
              fname = files(i).name;
            else
              fname = common_substring(fname, files(i).name);
            end
          end
        end
      end
    end

    parameters.ImageBaseName = fname;
  end

  if (isempty(parameters.ImageBaseName))
    if (nargout == 0)
      clear FileName;
    end

    return;
  end

  %Image format and number of characters for the images' number 
  ImagesList = dir(fullfile(parameters.ImageFolder, [parameters.ImageBaseName '*']));
  if isempty(ImagesList)
      error('In the selected ImageFolder there are not images with the selected ImageBaseName.')
  else
      StringNameImage1 = ImagesList(1,1).name;
      %PositionsUnderScores = strfind(StringNameImage1, '_');
      %PositionLastUnderScores = PositionsUnderScores(end);
      PositionEndCommon = length(parameters.ImageBaseName);
      PositionsPoints = strfind(StringNameImage1, '.');
      PositionLastPoint = PositionsPoints(end);
      ImageFormat = StringNameImage1(PositionLastPoint:end);
      parameters.ImageFormat = ImageFormat;
      %CharactersNumber = StringNameImage1(PositionLastUnderScores+1:PositionLastPoint-1);
      CharactersNumber = StringNameImage1(PositionEndCommon+1:PositionLastPoint-1);
      NumberCharactersNumber = ['%.' num2str(length(CharactersNumber)) 'd'];
      parameters.NumberCharactersNumber = NumberCharactersNumber;

      %clear StringNameImage1 PositionsUnderScores PositionLastUnderScores PositionsPoints
  end 
  if isempty(parameters.ImageIndexs)
      LengthImagesList = length(ImagesList);
      for i = 1:LengthImagesList
          StringNameImage1 = ImagesList(i,1).name;
          CurrentImageIndex = str2num(StringNameImage1(PositionLastPoint-length(CharactersNumber):PositionLastPoint-1));
          parameters.ImageIndexs = [parameters.ImageIndexs CurrentImageIndex];
      end
  end

  disp('------------------------------------------');
  disp('MicroMos: START.');

  % Determine the acquisition grid
  if (parameters.flag_GriddedAcquisition)
    disp('Determining the size of the acquisition grid...');
    [MatricesGLOBAL, mosaic_order] = DefineAcquisitionGrid(parameters);

    if (isempty(mosaic_order))
      disp('Troubles estimating the acqusition grid, trying with more corners...');

      tmp_params = parameters;
      tmp_params.numberCorners = 2*tmp_params.numberCorners;

      [MatricesGLOBAL, mosaic_order] = DefineAcquisitionGrid(tmp_params);
      if (isempty(mosaic_order))
        error('The acquisition grid could not be estimated using the provided data !');
      end
    end

    parameters.flag_ComputeRegistrations = 0;
    parameters.ImageIndexs = parameters.ImageIndexs(mosaic_order);
    MatricesGLOBAL = MatricesGLOBAL(:,:,mosaic_order);
  end

  start_index = 1;
  stop_index = length(parameters.ImageIndexs);

  %Vignetting function loading:
  Field = [];
  if parameters.flag_FlatField == 1
      disp('Computing the vignetting function...');
      LengthImagesList = length(ImagesList);
      for i = 1:LengthImagesList
         tmp_img = imread(fullfile(parameters.ImageFolder, ImagesList(i,1).name));
         Field = imvignette(tmp_img, Field);
      end
      Field = imvignette(Field);
      Field = fPixelAccuracy(Field);
      Field = Field./mean(Field(:));
      Field = Field(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
  end

  %% INITIALIZATION OF THE MOSAIC

  %Reference image loading and pre-processing
  strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(start_index));
  if strcmp(ImageFormat, '.mat')
      referenceFrame = load(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
      referenceFrame = cell2mat(struct2cell(referenceFrame));
  else
      referenceFrame = imread(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
  end
  referenceFrame = referenceFrame(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
  norm_factor = 1;
  if isa(referenceFrame, 'uint16'); norm_factor = 255/(2^16-1); end
  if isa(referenceFrame, 'uint32'); norm_factor = 255/(2^32-1); end
  if isa(referenceFrame, 'uint64'); norm_factor = 255/(2^64-1); end
  referenceFrame  = norm_factor*double(referenceFrame);
  Field = norm_factor*Field;

  if parameters.flag_Color==0
      if size(referenceFrame, 3)~=1
          referenceFrame = rgb2gray(referenceFrame);
      end
  else
      if size(referenceFrame, 3)==1
          parameters.flag_Color = 0;
          parameters.flag_LookUpTable = 0;
      end
  end
  referenceFrame = fPixelAccuracy(referenceFrame);
  if parameters.flag_FlatField == 1
      [referenceFrame, LookUpTable] = FlatFieldCorrection(referenceFrame, Field, parameters.flag_LookUpTable, LookUpTable);

      referenceFrame = fPixelAccuracy(referenceFrame);
  end

  %GLOBAL registration matrices loading:
  if parameters.flag_ComputeRegistrations ~= 1
      if (~parameters.flag_GriddedAcquisition)
        DirList = dir(['REGISTRATIONMATRIX' filesep '*.mat']);
        if isempty(DirList)
            error('In the folder called REGISTRATIONMATRIX there is not a file ".mat" related to the registration matrices.')
        else
            warning('If external matrices are loaded we suggest to do not use image rescaling to avoid errors.')
            Struttura = load(['REGISTRATIONMATRIX' filesep DirList(1).name]);
            copyfile(['REGISTRATIONMATRIX' filesep DirList(1).name], ['OUTPUT' filesep DirList(1).name]);
            RM = cell2mat(struct2cell(Struttura));
        end
        clear Struttura DirList
        MatricesGLOBAL = RM;
      end

      NumberOfregisteredImages = size(MatricesGLOBAL, 3);
      Indeces = [1:NumberOfregisteredImages];

      disp(['Registration of ' num2str(NumberOfregisteredImages) ' frames loaded.']);
  else

    %% MOSAIC REGISTRATION
    disp('Registering frames pairwise...')

    %Cycle for each image to be stitched
    base = referenceFrame;
    index = start_index;
    Indeces = index;
    NumberOfregisteredImages = 1;
    while index < stop_index

        if parameters.flag_SeekBestImages == 1
            % Automatic selection of the images to be registered.
            if index == start_index
                if parameters.flag_BleachingCorrection == 1
                    Vector_Indeces = [start_index+1, start_index+2];
                else
                    [Vector_RMSE, Vector_Indeces] = DefineBestImageToBeRegistered(parameters, index, stop_index, base, Field, LookUpTable, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                end
            else
                [Vector_RMSE, Vector_Indeces] = DefineBestImageToBeRegistered(parameters, index, stop_index, base, Field, LookUpTable, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
            end
            if size(Vector_Indeces)==1
                break
            end
            Indeces = [Indeces, Vector_Indeces(end-1)];
            TestNumber = 1;
        elseif parameters.flag_SeekBestImages == 0
            Vector_Indeces = index;
            Indeces = [Indeces, index+1];
            TestNumber = 0;
        end

        flag_Problem = 1; % if parameters.flag_SeekBestImages == 1 and a problem happens, the image to be registered is computed again when possible.
        while flag_Problem == 1

            if length(Vector_Indeces)<TestNumber+1
                % if a problem happened and parameters.flag_SeekBestImages == 0 or it is not possible to define another image to be registered, the mosaic updaiting stops.
                strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(end));
                disp(['STOP mosaic building: the algorithm is not able to find and image with a good overlap with: ' parameters.ImageBaseName strnum '.']);
                warning('STOP mosaic building: the algorithm is not able to find more images with a good overlap.')
                stop_index = index;
                Indeces = Indeces(1:end-1);
                ImageIndexs = parameters.ImageIndexs(Indeces);
                %save(['OUTPUT' filesep 'ImageIndexs.mat'], 'ImageIndexs');
                break
            end
            if parameters.flag_SeekBestImages == 1
                Indeces(end) = Vector_Indeces(end-TestNumber);
            end
            if index+1 == Indeces(end-1)
                % if a problem happened and parameters.flag_SeekBestImages == 0 or it is not possible to define another image to be registered, the mosaic updaiting stops.
                strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(end));
                disp(['STOP mosaic building: the algorithm is not able to find and image with a good overlap with: ' parameters.ImageBaseName strnum '.']);
                warning('STOP mosaic building: the algorithm is not able to find more images with a good overlap.')
                stop_index = index;
                Indeces = Indeces(1:end-1);
                ImageIndexs = parameters.ImageIndexs(Indeces);
                %save(['OUTPUT' filesep 'ImageIndexs.mat'], 'ImageIndexs');
                break
            end

            %Image to be stitched loading and pre-processing
            strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(Indeces(end)));
            if strcmp(parameters.ImageFormat, '.mat')
                unregistered = load(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
                unregistered = cell2mat(struct2cell(unregistered));
            else
                unregistered = imread(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
            end
            unregistered = unregistered(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
            unregistered = norm_factor*double(unregistered);
            if (parameters.flag_Color==0)
                if size(unregistered, 3)~=1
                    unregistered = rgb2gray(unregistered);
                end
            end
            unregistered = fPixelAccuracy(unregistered);
            if parameters.flag_FlatField == 1
                [unregistered, LookUpTable] = FlatFieldCorrection(unregistered, Field, parameters.flag_LookUpTable, LookUpTable);
                unregistered = fPixelAccuracy(unregistered);
            end

            %Registration estimation
            if parameters.flag_PhaseCorrelationOnly == 1
                if (parameters.flag_Color==0)
                    HF2F = RegistrationMatrixByPhaseCorrelationOnly(base, unregistered, parameters.PCscaleFactor, parameters.flag_PCglobalORlocal);
                else
                    ba = rgb2gray(uint8(base));
                    un = rgb2gray(uint8(unregistered));
                    ba = fPixelAccuracy(ba);
                    un = fPixelAccuracy(un);
                    HF2F = RegistrationMatrixByPhaseCorrelationOnly(ba, un, parameters.PCscaleFactor, parameters.flag_PCglobalORlocal);
                    clear ba un
                end
                disp(['Frame ' num2str(NumberOfregisteredImages) ': ' parameters.ImageBaseName strnum ' registered using the Phase Correlation Algorithm only. ']);
            else
                %% CORNER POINTS ESTIMATION
                numberCorners = 150;
                flag_Harris = 0;

                if (parameters.flag_Color==0)
                    [PointsBase, PointsTracked] = PointsDetectionTracking(base, unregistered, numberCorners, flag_Harris, parameters.ShiftEstimationMode, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                else
                    ba = rgb2gray(uint8(base));
                    un = rgb2gray(uint8(unregistered));
                    ba = fPixelAccuracy(ba);
                    un = fPixelAccuracy(un);
                    [PointsBase, PointsTracked] = PointsDetectionTracking(ba, un, numberCorners, flag_Harris, parameters.ShiftEstimationMode, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                    clear ba un
                end

                if (parameters.RegistrationMode == 0 && length(PointsBase) < 4) || (parameters.RegistrationMode == 1 && length(PointsBase) < 3) || (parameters.RegistrationMode == 2 && length(PointsBase) < 1)
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-frame registration: the overlapp between the image "' parameters.ImageBaseName strnum '" and the last one registered is too small.'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

                clear numberCorners flag_Harris

                %% FRAME-TO-FRAME REGISTRATION MATRIX ESTIMATION
                [HF2F, inliersF2F] = WarpingRegistrationMode(parameters.RegistrationMode, PointsBase, PointsTracked, parameters.RANSACerror);

                if isempty(HF2F) | size(HF2F) ~= [3, 3]
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-frame registration: the matrix estimated for the image '  parameters.ImageBaseName strnum ' is not valid.'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

                disp(['Frame ' num2str(NumberOfregisteredImages) ': ' parameters.ImageBaseName strnum '. ' 'Number tracked points/total points: ' num2str(length(inliersF2F)) '/' num2str(length(PointsBase))]);
                clear PointsBase PointsTracked
            end

            %% GLOBAL matrix updating
            GLOBAL = GLOBAL*inv(HF2F);

            %% FRAME WARPING AND STITCHING
            MatricesGLOBAL(:,:,end+1) = GLOBAL;

            % Next step
            base = unregistered;
            clear unregistered HF2F inliers2TF

            flag_Problem = 0;
        end
        index = Indeces(end);
        NumberOfregisteredImages = NumberOfregisteredImages + 1;
    end
  end

  [nrows, ncols, nchannels] = size(referenceFrame);
  is_rgb = (nchannels~=1);

  % MOSAIC CREATION
  disp('Building the mosaic...')

  NumberOfregisteredImages = length(Indeces);

  %Mosaic initialization to referenceFrame 
  [Mosaic, MosaicOrigin] = InitMosaic([nrows, ncols, nchannels], MatricesGLOBAL);

  for i=1:NumberOfregisteredImages

      % Get the registration matrix
      GLOBAL = MatricesGLOBAL(:,:,i);

      %Image to be stitched loading and pre-processing
      strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(Indeces(i)));
      disp(['Blending in ' parameters.ImageBaseName strnum '...']);

      if strcmp(parameters.ImageFormat, '.mat')
          unregistered = load(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
          unregistered = cell2mat(struct2cell(unregistered));
      else
          unregistered = imread(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum ImageFormat]));
      end
      unregistered = unregistered(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
      unregistered = norm_factor*double(unregistered);
      if (parameters.flag_Color==0)
          if is_rgb
              unregistered = rgb2gray(unregistered);
          end
      end
      unregistered = fPixelAccuracy(unregistered);
      if parameters.flag_FlatField == 1
          [unregistered, LookUpTable] = FlatFieldCorrection(unregistered, Field, parameters.flag_LookUpTable, LookUpTable);
          unregistered = fPixelAccuracy(unregistered);
      end

      %% BLEACHING CORRECTION
      %clear base
      %base = unregistered;
      if parameters.flag_BleachingCorrection == 1 && i>1
          for c = 1:nchannels
              [newunregistered, regionOverlapped] = ImagesForFTM(Mosaic(:,:,c), unregistered(:,:,c), GLOBAL, MosaicOrigin, parameters.RegistrationMode);
              LUT = BleachingLUTbuild(regionOverlapped, newunregistered);
              FrameOut = BleachingLUTuse(unregistered(:,:,c), LUT);
              unregistered(:,:,c) = FrameOut;
              clear regionOverlapped newunregistered LUT FrameOut
          end
      end

      %% FRAME-TO-MOSAIC REGISTRATION MATRIX ESTIMATION
      if parameters.flag_FrameToMosaic == 1 && i>1
          if (parameters.flag_Color==0)
              [newunregistered, regionOverlapped] = ImagesForFTM(Mosaic, unregistered, GLOBAL, MosaicOrigin, parameters.RegistrationMode);
          else
              newunregistered    = fPixelAccuracy(rgb2gray(uint8(unregistered)));
              regionOverlapped   = fPixelAccuracy(rgb2gray(uint8(Mosaic)));
              [newunregistered, regionOverlapped] = ImagesForFTM(regionOverlapped, newunregistered, GLOBAL, MosaicOrigin, parameters.RegistrationMode);
          end

          numberCorners = 150;
          flag_Harris = 0;

          [PointsBase, PointsTracked] = PointsDetectionTracking(double(regionOverlapped), double(newunregistered), numberCorners, flag_Harris, 0, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
          clear newunregistered regionOverlapped

          if (parameters.RegistrationMode == 0 && length(PointsBase) < 4) || (parameters.RegistrationMode == 1 && length(PointsBase) < 3) || (parameters.RegistrationMode == 2 && length(PointsBase) < 1)
              % a problem happened. If possible another image to be registered will be defined.
              disp(['Frame-to-mosaic registration: no enought matchings between the image "' parameters.ImageBaseName strnum '" and the mosaic.'])
          else

              [HF2M, inliersF2M] = WarpingRegistrationMode(parameters.RegistrationMode, PointsBase, PointsTracked, parameters.RANSACerror);
              clear PointsBase PointsTracked

              if isempty(HF2M) | size(HF2M) ~= [3, 3]
                  % a problem happened. If possible another image to be registered will be defined.
                  disp(['Frame-to-mosaic registration: the matrix estimated for the image '  parameters.ImageBaseName strnum ' is not valid.'])
              else
                  GLOBAL = GLOBAL*inv(HF2M);
                  MatricesGLOBAL(:,:,i) = GLOBAL;
              end
              clear HF2M inliersF2M
          end
      end

      %% Mosaic Updating
      [Mosaic] = MosaicUpdating(Mosaic, unregistered, GLOBAL, parameters.flag_Blending, MosaicOrigin, parameters.RegistrationMode);
  end

  goods = any(~isnan(Mosaic), 3);
  goodx = any(goods,1);
  goody = any(goods,2);

  Mosaic = Mosaic(goody, goodx, :);

  ImageIndexs = parameters.ImageIndexs(Indeces);
  parameters.ImageIndexs = ImageIndexs;
  parameters.Registrations = MatricesGLOBAL;
  %save(['OUTPUT' filesep 'ImageIndexs.mat'], 'ImageIndexs');

  %% MOSAIC POST PROCESSING: White Balancing.
  if parameters.flag_WhiteBalancing > 0
      if parameters.flag_Color == 1
          if parameters.flag_WhiteBalancing == 2          
              % Loading of the external RGB image to be used as colours reference.
              WB = [];
              DirList = dir(['WHITEBALANCING' filesep '*.*']);
              if isempty(DirList)
                  error('In the folder called WHITEBALANCING there is not a 3-channel image to be used as reference for the White Balancing.')
              else
                  WB = imread(['WHITEBALANCING' filesep DirList(end).name]);
                  copyfile(['WHITEBALANCING' filesep DirList(end).name], ['OUTPUT' filesep DirList(end).name]);
              end
              if size(WB,3)==3
                  % White Balancing of the output mosaic using the external RGB image as colours reference.
                  Mosaic = RGBWhiteBalancing(Mosaic, WB, -3);
              else
                  error('In the folder called WHITEBALANCING there is not a 3-channel image to be used as reference for the White Balancing.')
              end
          else
              % White Balancing of the output mosaic using the output mosaic itself as colours reference.
              Mosaic = RGBWhiteBalancing(Mosaic, Mosaic, -3);
          end
      end
  end

  % Convert back to the proper type, adjust and save
  if (norm_factor == 255/(2^16-1))
    Mosaic = uint16(Mosaic);
  elseif (norm_factor == 255/(2^32-1))
    Mosaic = uint32(Mosaic);
  elseif (norm_factor == 255/(2^64-1))
    Mosaic = uint64(Mosaic);
  else
    Mosaic = uint8(Mosaic);
  end

  if (parameters.flag_AdjustIntensityValues)
    Mosaic = imadjust(Mosaic, stretchlim(Mosaic));
  end

  fname = parameters.ImageFolder;
  if (fname(end) == filesep)
    fname = fname(1:end-1);
  end
  fname = [fname ImageFormat];

  imwrite(Mosaic, fname, ImageFormat(2:end));
  clear Mosaic

  if (nargout > 0)
    FileName{end+1} = fname;
  else
    clear FileName;
  end

  disp('MicroMos: THE END.');
end

function [Mosaic, MosaicOrigin] = InitMosaic(img_size, MatricesGLOBAL)
% Initializes the Mosaic, creating the array and setting up its origin

  nimgs = size(MatricesGLOBAL, 3);
  height = [0 0];
  width = [0 0];

  % Just go through all transformations and extract the outer most box
  for i=1:nimgs
    GLOBAL = MatricesGLOBAL(:,:,i);

    ULC=GLOBAL*[0;0;1];
    ULC=ULC./ULC(3);
    DLC=GLOBAL*[0;img_size(1)-1;1];
    DLC=DLC./DLC(3);
    DRC=GLOBAL*[img_size(2)-1;img_size(1)-1;1];
    DRC=DRC./DRC(3);
    URC=GLOBAL*[img_size(2)-1;0;1];
    URC=URC./URC(3);

    limits = [ULC DLC DRC URC];
    width(1) = min([limits(1,:), width(1)]);
    width(2) = max([limits(1,:), width(2)]);
    height(1) = min([limits(2,:), height(1)]);
    height(2) = max([limits(2,:), height(2)]);
  end

  % Add some spacing for safety
  width = width + nimgs*[-5 5];
  height = height + nimgs*[-5 5];

  % Create the arrays
  Mosaic = NaN([ceil(diff(height)) ceil(diff(width)) img_size(3)]);
  MosaicOrigin = -round([width(1) height(1)]);

  return;
end
