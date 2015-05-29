function parameters = ParametersByUser()
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ParametersByUser
% 
% Please read the file "_readme.txt", change the parameters into the file
% "ParametersByUser.mat" and then play the script "START.m" to build the
% mosaic. Please read carefully all the comments referring to the 
% parameters listed below.
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

%% PARAMETERS THAT MUST BE CHANGED BY THE USER
parameters.ImageFolder                      = '/Users/blanchou/Documents/MATLAB/Movies/Sectioning/Section5'; % absolute or relative path where the images are stored inside the computer. E.g.: 'INPUTIMAGES\' or 'C:\AbsolutePath\INPUTIMAGES\'.
parameters.ImageBaseName                    = 'Section5_'; % name of the images to be stitched, without the final cardinal number. The last character must be the underscore. E.g.: 'ImageMesenchymal_' for the image: "ImageMesenchymal_002.tif". The images must be saved (or copied) in the following format: "Name_#...##.format".
parameters.ImageIndexs                      = []; % arabic cardinal numbers of the specific images to be stitched between the ones present in the "ImageFolder". The numbers must be reported in the order of the images to be stitched. E.g.: [8, 72, 69] to process the images named "ImageMesenchymal_008.tif", "ImageMesenchymal_072.tif" and "ImageMesenchymal_069.tif" in this order.
parameters.flag_Color                       = 1; % (by default: 1). 0 to obtain the final mosaic in grey levels (also if the original images are RGB). 1 to obtain the final mosaic in RGB (only if the oiginal images are RGB). 
parameters.flag_SeekBestImages              = 0; % (by default: 1). The input frame are pre-selected to optimize the registration task. If this parameter is set to 1, not all the images pointed out are stitched in the final mosaic.
parameters.flag_BleachingCorrection         = 0; % (by default: 0 for bright field and phase contrast images, 1 for fluorescent images). To correct intensity decay, due to the photo-bleaching effect, between the images to be stitched. This phenomenon happen especially if the images to be stitched are fluorescent images. 0 = no bleaching correction; 1 = yes bleaching correction.
