function parameters = ParametersDefault(parameters)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ParametersDefault
% 
% Default parameters. We suggest to the users of do not change this code.
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

%% General parameters

%%% According to Piccinini, 2013, projective and affine give the same results, translative slightly worse BUT better if very little overlap
parameters.RegistrationMode                 = uint8(2); % (by default: uint8(2)). To choose the registration model that must be used to register the images. The registration model can be chosen between projective (suggested), affine and translative. Set: uint8(0) = 'projective' or uint8(1) = 'affine' or uint8(2) = 'translative' (computed at sub-pixel accuracy using the Lukas-Kanade feature tracker) or uint8(3) = 'translative' (computed at pixel accuracy using the phase-correlation algorithm only).

%%% Piccinini, 2013, Blending not required in widefield microscopy
parameters.flag_Blending                    = 0; % 0 = no blending; 1 = blending using a biquadratic function with the highest value in the image's centre (seams are only attenuated). 2 = linear blending with boundary values to completely avoid seams in the stitching regions (slow computation using pixel-based interpolation). 3 = linear blending with boundary values to completely avoid seams in the stitching regions (fast computation using Delaunay triangulation).

%%% No clue what this is.....
parameters.flag_WhiteBalancing              = 0; % 1 to perform the white balancing of the output mosaic using the mosaic itself as colour reference. 2 to perform the white balancing of the output mosaic loading an external 3-channel image (a RGB image) that must be copied in the folder called "WHITEBALANCING". 0 otherwise.

%%% Piccinini, 2012-2013, Crucial to get proper mosaic !!
parameters.flag_FlatField                   = 1; % 1 to flat-field correct the images using an input vignetting function. The vignetting function must be saved as matlab matrix in the folder named: "VIGNETTINGFUNCTION". In the "VIGNETTINGFUNCTION" folder must contain at maximum one vignetting function file.

%%% Piccinini, 2013, F2M clearly better !
parameters.flag_FrameToMosaic               = 1; % (by default: 1). 0 for registering the images according to the Frame-to-Frame registration approach; 1 (suggested) for registering the images according to the Frame-to-Mosaic registration approach.

%%% No clue what this is.....
parameters.RANSACerror                      = 2; % maximum subpixel reprojection error used inside RANSAC algorithm. We suggest 2 pixels.

% Registration mode selection
if parameters.RegistrationMode == 3;
    parameters.flag_PhaseCorrelationOnly    = 1; % It can assume values 0 or 1. 1 means that the images are registered according to the Phase Correlation ALgorithm only.
    parameters.RegistrationMode = 2;
    parameters.flag_FrameToMosaic = 0; % The registration mode is set to frame-to-frame when "parameters.RegistrationMode == uint8(3)".
else
    parameters.flag_PhaseCorrelationOnly = 0;
end

% Parameters used when the Phase Correlation is used to estimate the x-y translational shift between image pairs
parameters.flag_PCglobalORlocal             = 0; % It can assume values 0 or 1. 0 means that the metric used to determine the best shift inside the Phase Correlation ALgorithm is the global RMSE performed on the whole overlapping region. 1 means that the used metric is the RMSE performed using only the pixels with highest value. 
parameters.PCscaleFactor                    = 1; % to speed up the computational processes. Image rescale factor applyed to optionally resize the images. The value must be a positive integer. E.g.: ceil(2) to obtain half-sized images than the original images. 1 (suggested) means: rescaling not active.

% For the registration:
parameters.flag_ComputeRegistrations        = 1; % 0 for loading an external registration matrix to stitch the images. In case, the registration matrix must be saved as 3x3xn (n=number of images to be registered) Matlab matrix in the folder named: "REGISTRATIONMATRIX".
if parameters.flag_ComputeRegistrations == 0
    parameters.flag_SeekBestImages = 0;  % Do not change! There are not other possibilities implemented!
end

parameters.InterpolationMode                = 'bilinear'; % interpolation used to warp the image to be stitched. It can be: 'bicubic' or 'bilinear' (suggested) or 'nearest'.

% For correcting false-colours generated from the flat-field correction of colour images:
parameters.flag_LookUpTable = 1;
if parameters.flag_Color == 1
    if parameters.flag_FlatField == 1
        parameters.flag_LookUpTable = 1; % to use a 256 levels Look-Up-Table to map the grey levels into a single RGB conversion.
    end
end

% White Balancing
if parameters.flag_Color == 0
    parameters.flag_WhiteBalancing = 0;
end

% For partial-results visualization:
parameters.flag_ShowROIsingleImages         = 0; % to show inside the final mosaic the Region Of Interest of the single images. 0 = not active. 1 = active.
parameters.flag_SlideShow                   = 0; % to save an image of the current mosaic every time an image is stitched. 0 = not active. 1 = active.

% To assess the quality of the output mosaic using standard metrics computed only in the overlapping regions:
parameters.flag_ComputeMetrics              = 0; % to assess the quality of the output mosaic using standard metrics computed only in the overlapping regions. 0 = not active. 1 = active.

% To work with rescaled images:
parameters.ScaleFactor                      = ceil(1); % to speed up the computational processes. Image rescale factor applyed to optionally resize the images. The value must be a positive integer. E.g.: ceil(2) to obtain half-sized images than the original images. 1 (suggested) means: rescaling not active.
parameters.PixelAccuracy                    = 0; % (by default: 0). Pixel-accuracy of the output mosaic. 0 = double precision (slower computation, higher precision), 1 = single precision (faster computation, lower precision). If your computer have not enough memory to build the mosaic, we suggest to try to set this parameter at the value 1 (single precision).

%% CHECKS ON THE PARAMETERS: ABSOLUTELY DO NOT CHANGE FROM HERE:

% check on "flag_PhaseCorrelationOnly":
if parameters.flag_PhaseCorrelationOnly == 1
    parameters.flag_FrameToMosaic = 0; % Do not change! There are not other possibilities implemented!
    parameters.RegistrationMode = uint8(2);  % Do not change! There are not other possibilities implemented!
end

% check on "parameters.flag_WhiteBalancing": 
if parameters.flag_Color == 0
    parameters.flag_WhiteBalancing = 0;
end

% check on "parameters.flag_ComputeRegistrations":
if parameters.flag_ComputeRegistrations == 0
    parameters.RegistrationMode = uint8(0);
end

% check on "parameters.ScaleFactor":
if parameters.ScaleFactor == 0
    warning('The value of "parameters.ScaleFactor" was 0 and it has been changed to 1.')
    parameters.ScaleFactor = ceil(1);
end

% check on "parameters.ImageBaseName":
if ~strcmp(parameters.ImageBaseName(end), '_')
    parameters.ImageBaseName = strcat(parameters.ImageBaseName,'_');
end

% check on "parameters.ImageFolder":
%if isunix()
%    Slash = '/';
%else
%    Slash = '\';
%end;
Slash = filesep;
if ~strcmp(parameters.ImageFolder(end), Slash)
    parameters.ImageFolder = strcat(parameters.ImageFolder,Slash);
end

% OTHER INIZIALIZATIONS
rand('twister', 5490);
[status, message, messageid] = rmdir('OUTPUT','s');
mkdir('OUTPUT')
