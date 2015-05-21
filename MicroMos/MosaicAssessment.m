function  [MetricsTotal, MetricsSingleFrames] = MosaicAssessment(parameters, Mosaic, MaskOverlap, MatricesGLOBAL, MosaicOrigin, LookUpTable)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: MosaicAssessment
% 
% To assess the quality of the mosaic "Mosaic" using standard metrics (Mean
% Squared Error - MSE, Root Mean Squared Error - RMSE, Signal to Noise
% Ration - SNR, Universal Quality Index - UQI) computed only in the
% overlapping regions between the n images stitched into the mosaic
% according to the matrices saved in "MatricesGLOBAL".
%
% PARAMETERS:
%  parameters       Structure containing the parameters set into the files 
%                   "ParametersByUser.m" and "ParametersDefault.m".
%  Mosaic           Matrix of the current mosaic built stitching the images
%                   from 1 to n.
%  MaskOverlap      Matrix of the same size of "Mosaic". Each pixel reports
%                   the index of the last image that wrote in that specific 
%                   position.
%  MatricesGLOBAL   3x3xn-1 (n = number of images composing the mosaic) 
%                   registration matrices used to warp and stitch the last 
%                   n-1 images into the mosaic.
%  MosaicOrigin     x-y coordinate (x-y coordinate means the column-row
%                   coordinate) in the current "Mosaic" of the ULC (Up 
%                   Left Corner) of the first image stitched.
%  LookUpTable      256 levels Look-Up-Table to map the grey version of the 
%                   images to be stitched into a single RGB conversion.
%                   This Look-Up-Table is used only if 
%                   "parameters.flag_LookUpTable" = 1. 
%
% OUTPUT:
%  MetricsTotal     Structure containt the values of MSE, RMSE, SNR and UQI
%                   computed considering in the same time all the pixels of
%                   the mosaic where at least two input stitched images
%                   wrote.
%  MetricsSingleFrames  Structure containt n-1 values of MSE, RMSE, SNR and 
%                   UQI, computed reprojecting back the first n-1 images
%                   stitched into the mosaic and considering for each
%                   metrics evaluation only the pixels of the mosaic, inside
%                   tha bounding box of the reprojected image, where at
%                   least an image m with m>n wrote.
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

%% PARAMETERS SETTING

if nargin < 5
    parameters.flag_LookUpTable = 0;
end

fPixelAccuracy = @double;

[rowsMosaic, columnsMosaic, channelsMosaic] = size(Mosaic);
if channelsMosaic>1
    Mosaic = double(rgb2gray(uint8(Mosaic)));
end
Mosaic = feval(fPixelAccuracy, Mosaic);

if parameters.RegistrationMode == 2; % Translative
    RegistrationModel = 'affine';
elseif parameters.RegistrationMode == 1; %Sarebbe l'affine
    RegistrationModel = 'affine';
else
    RegistrationModel = 'projective';
end

if parameters.flag_FlatField == 1
    DirList = dir(['VIGNETTINGFUNCTION\' '*.mat']);
    if isempty(DirList)
        parameters.flag_FlatField = 0;
    else
        Struttura = load(['VIGNETTINGFUNCTION\' DirList(1).name]);
        copyfile(['VIGNETTINGFUNCTION\' DirList(1).name], ['OUTPUT\' DirList(1).name]);
        Field = cell2mat(struct2cell(Struttura));
        Field = feval(fPixelAccuracy, Field);
        Field = Field./mean(Field(:));
        Field = Field(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
    end
    clear Struttura DirList
end

start_index = 1;
stop_index = length(parameters.ImageIndexs);

PixelsMosaic = [];
PixelsUnregistered = [];

%% METRICS COMPUTATION
for index=1:stop_index-1

    %Image to be stitched loading and pre-processing
    strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(index));
    if strcmp(parameters.ImageFormat, '.mat')
        unregistered = load([parameters.ImageFolder parameters.ImageBaseName strnum parameters.ImageFormat]);
        unregistered = cell2mat(struct2cell(unregistered));
    else
        unregistered = imread([parameters.ImageFolder parameters.ImageBaseName strnum parameters.ImageFormat]);
    end
    unregistered = unregistered(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
    if isa(unregistered, 'uint16'); unregistered  = 255.*double(unregistered)./(2^16-1); end
    if isa(unregistered, 'uint32'); unregistered  = 255.*double(unregistered)./(2^32-1); end
    if isa(unregistered, 'uint64'); unregistered  = 255.*double(unregistered)./(2^64-1); end
    
    unregistered = feval(fPixelAccuracy, unregistered);
    
    if parameters.flag_FlatField == 1
        [unregistered LookUpTable] = FlatFieldCorrection(unregistered, Field, parameters.flag_LookUpTable, LookUpTable);
        unregistered = feval(fPixelAccuracy, unregistered);
    end
    if size(unregistered, 3)~=1
        unregistered = uint8(unregistered);
        unregistered = rgb2gray(unregistered);
        unregistered = feval(fPixelAccuracy, unregistered);
    end
    
    [rowsU, columnsU, channelsU] = size(unregistered);

    GLOBAL = MatricesGLOBAL(:,:,index);

    CurrentFrame = imtransform(unregistered, maketform(RegistrationModel,GLOBAL'), parameters.InterpolationMode, ...
        'XData',[MosaicOrigin(1) columnsMosaic-1+MosaicOrigin(1)],'YData',[MosaicOrigin(2) rowsMosaic-1+MosaicOrigin(2)],...
        'XYScale',[1],...
        'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
        'fill', NaN);

    [rowsPositions, colsPositions] = find(MaskOverlap>0 & MaskOverlap~=index-1 & isnan(CurrentFrame(:,:,1))~=1);

    if size(CurrentFrame,3)>1
        CurrentFrame = double(rgb2gray(uint8(CurrentFrame)));
    end

    if ~isempty(rowsPositions)

        CurMatMosaic = zeros(1, length(colsPositions));
        for px=1:length(colsPositions)
             CurMatMosaic(px) = Mosaic(rowsPositions(px),colsPositions(px));
        end
        PixelsMosaic = [PixelsMosaic CurMatMosaic];

        CurUnreg = zeros(1, length(rowsPositions));
        for px=1:length(colsPositions)
             CurUnreg(px) = CurrentFrame(rowsPositions(px),colsPositions(px));
        end
        PixelsUnregistered = [PixelsUnregistered CurUnreg];

        MetricsSingleFrames{index} = Metrics(CurMatMosaic, CurUnreg);

        clear CurMatMosaic CurUnreg
    else
        MetricsSingleFrames{index} = [];
    end
    clear unregistered CurrentFrame
end

MetricsTotal = Metrics(PixelsMosaic, PixelsUnregistered);



function ms = Metrics(FrameOverlap, MosaicOverlap, Mask)
ms.MSE =  Metric_MSE(FrameOverlap, MosaicOverlap);
ms.RMSE = Metric_RMSE(FrameOverlap, MosaicOverlap);
ms.SNR =  Metric_SNR(FrameOverlap, MosaicOverlap);
ms.UQI =  Metric_UQI(FrameOverlap, MosaicOverlap);