function [unregistered, Mosaic] = ImagesForFTM(Mosaic, unregistered, GLOBAL, MosaicOrigin, RegistrationMode)
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

% Warping
unregistered = myimtransform(unregistered, modello, GLOBAL, [columnsMosaic rowsMosaic], MosaicOrigin);

% Current location of the warped image
warped = ~isnan(unregistered);
warpedx = any(warped,1);
warpedy = any(warped,2);

% Sub-image
unregistered     = unregistered(warpedy,warpedx,:);
Mosaic    = Mosaic(warpedy,warpedx,:);
