function  [OverlapPercentage, RegistrationErrors] = RegistrationErrorAssessment(base, unregistered, GLOBAL, RegistrationMode, InterpolationMode)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: RegistrationErrorAssessment
% 
% To assess the quality of the registrastration stage computing standard
% metrics (Mean Squared Error - MSE, Root Mean Squared Error - RMSE, Signal
% to Noise Ration - SNR, Universal Quality Index - UQI) only in the
% overlapping regions between the "unregistered" image to be stitched into
% "base" according to the matrices saved in "GLOBAL".
%
% PARAMETERS:
%  base             Matrix of the first image.
%  unregistered     Image to be stitched into "base".
%  GLOBAL           3x3x1 registration matrices used for stitching 
%                   "unregistered" into the base.
%  RegistrationMode     Modality used to warp the image "unregistered" to be 
%                   stitched into "base". uint8(0) = 'projective' or 
%                   uint8(1) = 'affine' or uint8(2) = 'translative'.
%  InterpolationMode    interpolation used to warp the image to be 
%                   stitched. It can be: 'bicubic' or 'bilinear' or 
%                   'nearest'.
%
% OUTPUT:
%  OverlapPercentage     Percentage of overlap between "unregistered" and
%                   "base"
%  RegistrationErrors  Structure containt the values of MSE, RMSE, SNR and 
%                   UQI, computed stitching the image "unregistered" to be
%                   stitched into "base".
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

[rowsbase, columnsbase, channelsbase] = size(base);
[rowsU, columnsU, channelsU] = size(unregistered);
if channelsbase ~= channelsU
    error('Input images of different size')
end
if channelsbase>1 && channelsU>1
    base = double(rgb2gray(uint8(base)));
    unregistered = double(rgb2gray(uint8(unregistered)));
end

if RegistrationMode == 2; % Translative
    RegistrationModel = 'affine';
elseif RegistrationMode == 1; %Sarebbe l'affine
    RegistrationModel = 'affine';
else
    RegistrationModel = 'projective';
end

baseOrigin = [0, 0];

CurrentFrame = imtransform(unregistered, maketform(RegistrationModel,GLOBAL'), InterpolationMode, ...
    'XData',[baseOrigin(1) columnsbase-1+baseOrigin(1)],'YData',[baseOrigin(2) rowsbase-1+baseOrigin(2)],...
    'XYScale',[1],...
    'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
    'fill', NaN);

%[rowsPositions, colsPositions] = find(CurrentFrame(:,:,1)>=0);
Positions = find(CurrentFrame(:,:,1)>=0);
OverlapPercentage = length(Positions)/(rowsbase*columnsbase);
CurMatbase = base(Positions);
CurUnreg = CurrentFrame(Positions);
RegistrationErrors = Metrics(CurMatbase, CurUnreg);

function ms = Metrics(FrameOverlap, MosaicOverlap, Mask)
ms.MSE =  Metric_MSE(FrameOverlap, MosaicOverlap);
ms.RMSE = Metric_RMSE(FrameOverlap, MosaicOverlap);
ms.SNR =  Metric_SNR(FrameOverlap, MosaicOverlap);
ms.UQI =  Metric_UQI(FrameOverlap, MosaicOverlap);