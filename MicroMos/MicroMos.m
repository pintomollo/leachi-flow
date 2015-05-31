function [Mosaic, MaskOverlap, MatricesGLOBAL] = MicroMos(parameters)
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

if (nargin == 0 || isempty(parameters))
  parameters = get_struct('MicroMos');
end

%% PARAMETERS SETTTING
GLOBAL = eye(3,3);
MatricesGLOBAL = GLOBAL;
MaskOverlap = [];
MosaicOrigin = [0 0];
LookUpTable = NaN.*ones(4,256);

if parameters.PixelAccuracy == 0
    fPixelAccuracy = @double;
else
    fPixelAccuracy = @single;
end

if (isempty(parameters.ImageFolder) || isempty(parameters.ImageBaseName))
  parameters.ImageFolder = uigetdir('Movies','Select the directory containing the mosaic');
  files = dir(fullfile(parameters.ImageFolder, '*'));

  fname = '';
  for i=1:length(files)
    if (files(i).name(1)~='.')
      if (isempty(fname))
        fname = files(i).name;
      else
        fname = common_substring(fname, files(i).name);
      end
    end
  end

  parameters.ImageBaseName = fname;
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

    %{
    DirList = dir(['VIGNETTINGFUNCTION' filesep '*.mat']);
    if isempty(DirList)
        disp('In the folder called VIGNETTINGFUNCTION there is not a file ".mat" related to the vignetting function.')
        disp('The input parameter "flag_FlatField" has been changed to 0.')
        parameters.flag_FlatField = 0;
    else
        Struttura = load(['VIGNETTINGFUNCTION' filesep DirList(1).name]);
        copyfile(['VIGNETTINGFUNCTION' filesep DirList(1).name], ['OUTPUT' filesep DirList(1).name]);
        Field = cell2mat(struct2cell(Struttura));
        Field = fPixelAccuracy(Field);
        Field = Field./mean(Field(:));
        Field = Field(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
    end
    clear Struttura DirList
    %}
end

%GLOBAL registration matrices loading:
if parameters.flag_ComputeRegistrations ~= 1
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

%{
if isa(referenceFrame, 'uint16'); referenceFrame  = 255.*double(referenceFrame)./(2^16-1); end
if isa(referenceFrame, 'uint32'); referenceFrame  = 255.*double(referenceFrame)./(2^32-1); end
if isa(referenceFrame, 'uint64'); referenceFrame  = 255.*double(referenceFrame)./(2^64-1); end
%}
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
    [referenceFrame LookUpTable] = FlatFieldCorrection(referenceFrame, Field, parameters.flag_LookUpTable, LookUpTable);
    referenceFrame = fPixelAccuracy(referenceFrame);
end

%Mosaic initialization to referenceFrame 
Mosaic = referenceFrame;
Corner_Position = [[0,0,1]',[size(Mosaic,2)-1,0,1]',[size(Mosaic,2)-1,size(Mosaic,1)-1,1]',[0,size(Mosaic,1)-1,1]',[MosaicOrigin(1), MosaicOrigin(2),1]'];

%% MOSAIC UPDATING 

%Cycle for each image to be stitched
base = referenceFrame;
clear referenceFrame
index = start_index;
Indeces = index;
disp('MicroMos: START.');
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
        %{
        if isa(unregistered, 'uint16'); unregistered  = 255.*double(unregistered)./(2^16-1); end
        if isa(unregistered, 'uint32'); unregistered  = 255.*double(unregistered)./(2^32-1); end
        if isa(unregistered, 'uint64'); unregistered  = 255.*double(unregistered)./(2^64-1); end
        %}
        if (parameters.flag_Color==0)
            if size(unregistered, 3)~=1
                unregistered = rgb2gray(unregistered);
            end
        end
        unregistered = fPixelAccuracy(unregistered);
        if parameters.flag_FlatField == 1
            [unregistered LookUpTable] = FlatFieldCorrection(unregistered, Field, parameters.flag_LookUpTable, LookUpTable);
            unregistered = fPixelAccuracy(unregistered);
        end

        %Registration estimation or loading
        if parameters.flag_ComputeRegistrations ~= 1
            GLOBAL = RM(:,:,index-start_index+1 +1);
            disp(['Frame ' num2str(NumberOfregisteredImages) ': ' parameters.ImageBaseName strnum ' registered.']);
        else
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
                flag_PhaseCorrelation = 1;

                if (parameters.flag_Color==0)
                    [PointsBase, PointsTracked] = PointsDetectionTracking(base, unregistered, numberCorners, flag_Harris, flag_PhaseCorrelation, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                else
                    ba = rgb2gray(uint8(base));
                    un = rgb2gray(uint8(unregistered));
                    ba = fPixelAccuracy(ba);
                    un = fPixelAccuracy(un);
                    [PointsBase, PointsTracked] = PointsDetectionTracking(ba, un, numberCorners, flag_Harris, flag_PhaseCorrelation, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                    clear ba un
                end

                if (parameters.RegistrationMode == 0 && length(PointsBase) < 4) || (parameters.RegistrationMode == 1 && length(PointsBase) < 3) || (parameters.RegistrationMode == 2 && length(PointsBase) < 1)
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-frame registration: the overlapp between the image "' parameters.ImageBaseName strnum '" and the last one registered is too small.'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

                clear numberCorners flag_Harris flag_PhaseCorrelation decFactorPC

                %% FRAME-TO-FRAME REGISTRATION MATRIX ESTIMATION
                try        
                    [HF2F, inliersF2F] = WarpingRegistrationMode(parameters.RegistrationMode, PointsBase, PointsTracked, parameters.RANSACerror);
                catch ME1
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-frame registration: problem in the model parameters estimation. It tryed to register the image: "' parameters.ImageBaseName strnum '" .'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

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

            %[OverlapPercentage(NumberOfregisteredImages, 1), RegistrationErrors{NumberOfregisteredImages}] = RegistrationErrorAssessment(base, unregistered, inv(HF2F), parameters.RegistrationMode, parameters.InterpolationMode);
            
            %% BLEACHING CORRECTION
            unregisteredOld = unregistered;
            if parameters.flag_BleachingCorrection == 1
                for c = 1:size(Mosaic,3)
                    [newunregistered, regionOverlapped] = ImagesForFTM(base(:,:,c), unregistered(:,:,c), inv(HF2F), [0, 0], parameters.InterpolationMode, parameters.RegistrationMode);        
                    LUT = BleachingLUTbuild(regionOverlapped, newunregistered);
                    FrameOut = BleachingLUTuse(unregistered(:,:,c), LUT);
                    unregistered(:,:,c) = FrameOut;
                    clear regionOverlapped newunregistered LUT FrameOut
                end
            end
            %clear base
            base = unregisteredOld;
            clear unregisteredOld
            
            
            %% GLOBAL matrix updating
            GLOBAL = GLOBAL*inv(HF2F);
            clear HF2F inliers2TF

            %% FRAME-TO-MOSAIC REGISTRATION MATRIX ESTIMATION
            if parameters.flag_FrameToMosaic == 1
                if (parameters.flag_Color==0)
                    [newunregistered, regionOverlapped] = ImagesForFTM(Mosaic, unregistered, GLOBAL, MosaicOrigin, parameters.InterpolationMode, parameters.RegistrationMode);        
                else
                    unregisteredGrey    = rgb2gray(uint8(unregistered));      
                    MosaicGrey          = rgb2gray(uint8(Mosaic));
                    unregisteredGrey    = fPixelAccuracy(unregisteredGrey);
                    MosaicGrey          = fPixelAccuracy(MosaicGrey);
                    [newunregistered, regionOverlapped] = ImagesForFTM(MosaicGrey, unregisteredGrey, GLOBAL, MosaicOrigin, parameters.InterpolationMode, parameters.RegistrationMode);
                    clear unregisteredGrey MosaicGrey
                end

                numberCorners = 150;
                flag_Harris = 0;
                flag_PhaseCorrelation = 0;

                [PointsBase, PointsTracked] = PointsDetectionTracking(double(newunregistered), double(regionOverlapped), numberCorners, flag_Harris, flag_PhaseCorrelation, parameters.flag_PCglobalORlocal, parameters.PCscaleFactor);
                clear newunregistered regionOverlapped
                clear numberCorners flag_Harris flag_PhaseCorrelation decFactorPC

                if (parameters.RegistrationMode == 0 && length(PointsBase) < 4) || (parameters.RegistrationMode == 1 && length(PointsBase) < 3) || (parameters.RegistrationMode == 2 && length(PointsBase) < 1)
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-mosaic registration: no enought matchings between the image "' parameters.ImageBaseName strnum '" and the mosaic.'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

                try        
                    [HF2M, inliersF2M] = WarpingRegistrationMode(parameters.RegistrationMode, PointsBase, PointsTracked, parameters.RANSACerror);
                catch ME1
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-mosaic registration: problem in the model parameters estimation. It tryed to register the image: "' parameters.ImageBaseName strnum '" .'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end
                clear PointsBase PointsTracked

                if isempty(HF2M) | size(HF2M) ~= [3, 3]
                    % a problem happened. If possible another image to be registered will be defined.
                    disp(['Frame-to-mosaic registration: the matrix estimated for the image '  parameters.ImageBaseName strnum ' is not valid.'])
                    flag_Problem = 1;
                    TestNumber = TestNumber + 1;
                    continue
                end

                GLOBAL = GLOBAL*HF2M;
                clear HF2M inliersF2M
            end
        end

        %% FRAME WARPING AND STITCHING
        MatricesGLOBAL(:,:,end+1) = GLOBAL;
        
        %% BLEACHING CORRECTION
        %clear base
        %base = unregistered;
        if parameters.flag_BleachingCorrection == 1
            for c = 1:size(Mosaic,3)
                [newunregistered, regionOverlapped] = ImagesForFTM(Mosaic(:,:,c), unregistered(:,:,c), GLOBAL, MosaicOrigin, parameters.InterpolationMode, parameters.RegistrationMode);        
                LUT = BleachingLUTbuild(regionOverlapped, newunregistered);
                FrameOut = BleachingLUTuse(unregistered(:,:,c), LUT);
                unregistered(:,:,c) = FrameOut;
                clear regionOverlapped newunregistered LUT FrameOut
            end
        end

        %% Mosaic Updating 
        
        % Memory check
        %{
        [user sys] = memory;
        CurrentMemUsedMATLAB = user.MemUsedMATLAB;
        CurrentMaxPossibleArrayBytes = user.MaxPossibleArrayBytes;
        CurrentMemAvailableAllArrays = user.MemAvailableAllArrays;
        MemoryRequiredForCopyesOfMosaic = 2; %Number of copyes of the matrix called Mosaic
        MosaicParameters = whos('Mosaic');
        MosaicMemory = MosaicParameters.bytes;
        if MosaicMemory < CurrentMaxPossibleArrayBytes && CurrentMemAvailableAllArrays >= (CurrentMemUsedMATLAB + (MemoryRequiredForCopyesOfMosaic*MosaicMemory))        
        %}
            % The most computational expensive function is the following one. Optimization would be necessary.
            [Mosaic, MosaicOrigin, MaskOverlap, Corner_Position] = MosaicUpdating(Mosaic, unregistered, GLOBAL, parameters.flag_Blending, MosaicOrigin, MaskOverlap, parameters.InterpolationMode, parameters.RegistrationMode, NumberOfregisteredImages, Corner_Position);
        %{
            clear MosaicMemory CurrentMemUsedMATLAB
        else
            disp(['STOP mosaic building: MemAvailableAllArrays is not sufficient for registering the image: ' parameters.ImageBaseName strnum]);
            warning('STOP mosaic building: MemAvailableAllArrays is not sufficient. Try to defragment your hard disk.')
            stop_index = index;
            ImageIndexs = parameters.ImageIndexs(Indeces);
            %parameters.ImageIndexs = ImageIndexs;
            save('OUTPUT\ImageIndexs.mat', 'ImageIndexs');
            break
        end
        %}

        %Show the final mosaic
        %figure(1), imshow(uint8(Mosaic), 'Border', 'Tight')

        clear unregistered        
        flag_Problem = 0;
    end
    index = Indeces(end);
    NumberOfregisteredImages = NumberOfregisteredImages + 1;
end

ImageIndexs = parameters.ImageIndexs(Indeces);
parameters.ImageIndexs = ImageIndexs;
%save(['OUTPUT' filesep 'ImageIndexs.mat'], 'ImageIndexs');
    
if parameters.flag_ShowROIsingleImages == 1 || parameters.flag_SlideShow == 1
%     Corner_Position = n times [Up Left Corner, U. Right C., Low R. C.,
%     LLC, "MosaicOrigin"]. n = number of images composing the mosaic.
%     Inside "Corner_Position" are saved the coordinated of the [ULC URC
%     LRC LLC "MosaicOrigin"] of the single n images stitched into the
%     mosaic. Here the coordinates are correctly scaled (+1 pixels) to find
%     the corner coordinates inside the matrix of the mosaic. 
    CornerPositionScaled = Corner_Position - [Corner_Position(1,end).*ones(1,size(Corner_Position, 2)); Corner_Position(2,end).*ones(1,size(Corner_Position, 2)); 0.*ones(1,size(Corner_Position, 2))] + [ones(1,size(Corner_Position, 2)); ones(1,size(Corner_Position, 2)); 0.*ones(1,size(Corner_Position, 2))];
    clear Corner_Position
end

%To show inside the final mosaic the Region Of Interest of the single images. 0 = not active.
if parameters.flag_ShowROIsingleImages == 1
    MosaicStretched = uint8(Mosaic); 
    if size(Mosaic, 3)== 1; pos = find(isnan(Mosaic)==0); MosaicStretched = uint8(imadjust(MosaicStretched,stretchlim(Mosaic(pos)),[])); end
    figure, imshow(MosaicStretched,'Border','Tight')
    hold on
    LengthCornerPositionScaled = size(CornerPositionScaled, 2);
    for i = 1:5:LengthCornerPositionScaled ;    
        hsvcolor=[0.6 0.6 0.6]+0.4*rand(1,3);
        colorline=[hsvcolor(1) hsvcolor(2) hsvcolor(3)];
        plot([CornerPositionScaled(1,i) CornerPositionScaled(1,i+1)], [CornerPositionScaled(2,i) CornerPositionScaled(2,i+1)],'LineWidth',2,'Color',colorline);
        plot([CornerPositionScaled(1,i+1) CornerPositionScaled(1,i+2)], [CornerPositionScaled(2,i+1) CornerPositionScaled(2,i+2)],'LineWidth',2, 'Color',colorline)
        plot([CornerPositionScaled(1,i+2) CornerPositionScaled(1,i+3)], [CornerPositionScaled(2,i+2) CornerPositionScaled(2,i+3)],'LineWidth',2, 'Color',colorline)
        plot([CornerPositionScaled(1,i+3) CornerPositionScaled(1,i)], [CornerPositionScaled(2,i+3) CornerPositionScaled(2,i)],'LineWidth',2, 'Color',colorline)
        clear hsvcolor colorline
    end
    hold off
    clear LengthCornerPositionScaled MosaicStretched
end

%To save an image of the current mosaic every time an image is stitched
if parameters.flag_SlideShow == 1    
    SlideShow(parameters, Mosaic, MatricesGLOBAL, MosaicOrigin, LookUpTable)
end

%% MOSAIC EVALUATION
if parameters.flag_ComputeMetrics == 1
    if parameters.flag_ComputeRegistrations == 1
        for i = 1:length(OverlapPercentage); OverlapRMSE(1:2,i) = [OverlapPercentage(i) RegistrationErrors{i}.RMSE]; end;
        OverlapRMSEmean = mean(OverlapRMSE');
        OverlapRMSEstd = std(OverlapRMSE');
        %save(['OUTPUT' filesep 'REGISTRATIONERROR.mat'], 'OverlapPercentage', 'RegistrationErrors', 'OverlapRMSE', 'OverlapRMSEmean', 'OverlapRMSEstd');
    end
    [MetricsTotal, MetricsSingleFrames] = MosaicAssessment(parameters, Mosaic, MaskOverlap, MatricesGLOBAL, MosaicOrigin, LookUpTable);
    numOfFrame = length(MetricsSingleFrames);
    for i = 1:numOfFrame
        ValuesTable(i,1) = MetricsSingleFrames{i}.MSE;
        ValuesTable(i,2) = MetricsSingleFrames{i}.RMSE;
        ValuesTable(i,3) = MetricsSingleFrames{i}.SNR;
        ValuesTable(i,4) = MetricsSingleFrames{i}.UQI;
    end
    ValuesMean  = mean(ValuesTable);
    ValuesStd   = std(ValuesTable);
    %save(['OUTPUT' filesep 'METRICS.mat'], 'MetricsTotal', 'MetricsSingleFrames', 'ValuesTable', 'ValuesMean', 'ValuesStd');
end

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

disp('MicroMos: THE END.');
