function [newunregistered, regionOverlapped] = ImagesForFTM(Mosaic, unregistered, GLOBAL, MosaicOrigin, InterpolationMode, RegistrationMode)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ImagesForFTM
% 
% To define the overlapping and the overlapped part when the image n
% "unregistered" is stitched on the current mosaic "Mosaic" built with the
% images from 1 to n-1. "GLOBAL" is the 3x3 transformation matrix to warp
% and stitch the image "unregistered" on the mosaic "Mosaic".
% "MosaicOrigin" is the x-y coordinate (x-y coordinate means the column-row
% coordinate) in the current "Mosaic" of the Up Left Corner of the first
% image stitched into the mosaic. "InterpolationMode" can be 'nearest' or
% 'bilinear' or 'bicubic' and point out the registration modality used to
% warp the images. "RegistrationMode" can be 0 (projective) or 1 (affine)
% or 2 (translative).
%
% PARAMETERS:
%  Mosaic           Matrix of the current mosaic built stitching the images
%                   from 1 to n-1.
%  unregistered     Image n to be stitched in the current mosaic.
%  GLABAL           3x3 registration matrix to warp and stitch the image
%                   "unregisterd" into the mosaic "Mosaic".
%  MosaicOrigin     x-y coordinate (x-y coordinate means the column-row
%                   coordinate) in the current "Mosaic" of the ULC (Up 
%                   Left Corner) of the first image stitched.
%  InterpolationMode    Interpolation used to warp the image to be 
%                   stitched. It can be: 'bicubic' or 'bilinear'
%                   (suggested) or 'nearest'.
%  RegistrationMode     To pre-fix the registration model that must be 
%                   used to register the image "unregistered". The 
%                   registration model can be chosen between 0 (projective,
%                   suggested) or 1 (affine) or 2 (translative).
%
% OUTPUT:
%  newunregistered  Version of the image "unregistered" warped according to
%                   the matrix "GLOBAL". The image is inscribed into the
%                   maximum non-NaN bounding box.
%  regionOverlapped Region of the current mosaic "Mosaic" that will be
%                   overlapped by "newunregistered". "regionOverlapped" and
%                   "newunregistered" are matrices of the same size.
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

[rowsMosaic, columnsMosaic, channelsMosaic] = size(Mosaic);
[rowsU, columnsU, channelsU] = size(unregistered);

if RegistrationMode == 2; % Translative
    modello = 'affine';
elseif RegistrationMode == 1;
    modello = 'affine';
else
    modello = 'projective';
end

% U = Up; D = Down; L = Left; R = Right; C = Corner. 
ULC=GLOBAL*[0;0;1];
ULC=ULC./ULC(3);
DLC=GLOBAL*[0;rowsU-1;1];
DLC=DLC./DLC(3);
DRC=GLOBAL*[columnsU-1;rowsU-1;1];
DRC=DRC./DRC(3);
URC=GLOBAL*[columnsU-1;0;1];
URC=URC./URC(3);

%Maximum Bounding Box
Xmin = (min([MosaicOrigin(1),ULC(1),DLC(1)])); %ceil
Xmax = (max([columnsMosaic-1+MosaicOrigin(1),URC(1),DRC(1)])); % floor
Ymin = (min([MosaicOrigin(2),ULC(2),URC(2)])); % ceil
Ymax = (max([rowsMosaic-1+MosaicOrigin(2),DLC(2),DRC(2)])); % floor
XminI = floor(Xmin);
XmaxI = ceil(Xmax);
YminI = floor(Ymin);
YmaxI = ceil(Ymax);

unregisteredWarped = imtransform(double(unregistered), maketform(modello,GLOBAL'), InterpolationMode, ...
    'XData',[XminI XmaxI],'YData',[YminI YmaxI],...
    'XYScale',[1],...
    'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
    'fill', NaN);    
clear unregistered

dx = XminI - MosaicOrigin(1);
dy = YminI - MosaicOrigin(2);
ox = dx;
oy = dy;

%Mosaic copied in the new Bounding Box. Always "nearest" interpolation.
newMosaic = imtransform(Mosaic, maketform('affine',eye(3)), 'nearest',...
    'XData',[ox XmaxI-XminI+ox],'YData',[oy YmaxI-YminI+oy],...
    'XYScale',[1],...
    'UData',[0 columnsMosaic-1],'VData',[0 rowsMosaic-1],...
    'fill', NaN); 

[ro, co] = find(isnan(unregisteredWarped)==0);
ULCrc = [min(ro), min(co)];
DRCrc = [max(ro), max(co)];

newunregistered     = unregisteredWarped(ULCrc(1):DRCrc(1), ULCrc(2):DRCrc(2), :);
regionOverlapped    = newMosaic(ULCrc(1):DRCrc(1), ULCrc(2):DRCrc(2), :);