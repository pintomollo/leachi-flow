% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 June 2013
% NAME: START
% 
% Run this script to obtain a mosaic starting from an image video.
% The input video must be recorded with no compression and no interlacing.
%
% Please read the file "_readme.txt", change the parameters into the file
% "ParametersByUser.mat" and then play this script "START.m" to build the
% mosaic.
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

%% WORKSPACE CLEANING
clc, clear, close all

%% PARAMETERS SET UP: the user must change this code below.
parameters = ParametersByUser();

%% OTHER PARAMETERS: We suggest to the user of do not change this code.
parameters = ParametersDefault(parameters);

%test_vignette(parameters);

%% MOSAIC BUILDING
[Mosaic, MaskOverlap, MatricesGLOBAL] = MicroMos(parameters);

%% OUTPUT SAVING
copyfile('ParametersByUser.m', ['OUTPUT' filesep 'ParametersByUser.m']);
copyfile('ParametersDefault.m', ['OUTPUT' filesep 'ParametersDefault.m']);
save(['OUTPUT' filesep 'parameters.mat'], 'parameters');
save(['OUTPUT' filesep 'Mosaic.mat'], 'Mosaic');
save(['OUTPUT' filesep 'MaskOverlap.mat'], 'MaskOverlap');
save(['OUTPUT' filesep 'MatricesGLOBAL.mat'], 'MatricesGLOBAL');
if size(Mosaic, 3)== 1; pos = find(isnan(Mosaic)==0); Mosaic = uint8(Mosaic); Mosaic = uint8(imadjust(Mosaic,stretchlim(Mosaic(pos)),[])); clear pos; end
imwrite(uint8(Mosaic), ['OUTPUT' filesep 'Mosaic.tif'], 'tif');
figure(), imshow(uint8(Mosaic), 'Border', 'Tight');
