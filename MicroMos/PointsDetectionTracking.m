function [PointsBase, PointsTracked] = PointsDetectionTracking(base, unregistered, numberCorners, flag_Harris, flag_PhaseCorrelation, flag_PCglobalORlocal, PCscaleFactor) 
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: PointsDetectionTracking
% 
% To estimate correspondences between the 2D matrices "base" and
% "unregistered". "base" and "unregistered" must be matrices of the same
% dimension. Two different types of features can be extracted. If
% "flag_Harris" is 1 the Harris points are extracted, otherwise the
% Shi-Tomasi corner points. In particular, "numberCorners" points are
% extracted from "base" and tracked in "unregistered". If
% "flag_PhaseCorrelation" is 1 the Phase Correlation algorithm is use to
% estimate a guess shift between "base" and "unregistered". If
% "PCscaleFactor" is a positive integer number greater than 1, the Phase
% Correlation is computed using a rescaled version of "base" and
% "unregistered". "PointsBase" and "PointsTracked" are the x-y coordinates
% (x-y means the column-row coordinate inside the 2D matrix) of the points
% extracted in "base" and tracked in "unregistered", respectively.
%
% PARAMETERS:
%  base             Input image used to extract the corner points. We
%                   suggest to cast to double and convert to mono-channel
%                   image.
%  unregistered     Input image used to track the corner points extracted
%                   from "base". We suggest to cast to double and convert 
%                   to mono-channel image.
%  numberCorners    maximum number of corners extracted from "base". We   
%                   suggest 150. 
%  flag_Harris      1 to extract Harris corner points. 0 to extract
%                   Shy-tomasi corner points. We suggest 0.
%  flag_PhaseCorrelation    1 to use the Phase Correlation algorithm to 
%                   compute an initial guess for the tracking of the
%                   corners.
%  flag_PCglobalORlocal  It can assume values 0 or 1. 0 means that the 
%                   metric used to determine the best shift inside the 
%                   Phase Correlation ALgorithm is the global RMSE 
%                   performed on the entire overlapping regions. 1 means 
%                   that the used metric is the RMSE performed using only 
%                   the pixels with highest value. 
%  PCscaleFactor    The value must be a positive integer. If greater than 1
%                   the input images are rescaled before computing the 
%                   Phase Correlation.
%
% OUTPUT:
%  PointsBase       2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  PointsTracked    2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
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

% The corner points are extracted in the "base" image and they are tracked
% in the "unregistered" image. The functions ("STFeatureExtract.mexw32" and
% "mexFunctionLKT.mexw32") used to registered the images have been built by
% importing functions of the C++ OpenCV (see file "_licenseOpenCV.txt") and
% they require that the files with extension ".dll" are already present in
% the main folder.

if ~all(size(base)==size(unregistered))
    error('The two input matrices must be of the same dimension.');
end

[rows, columns, channels] = size(base);

%% Features extraction

if (flag_Harris)
  method = 'Harris';
else
  method = 'MinimumEigenvalue';
end
tmp_base = double(base);
tmp_base(isnan(base)) = nanmean(tmp_base(:));
detection_points2 = corner(tmp_base, method, numberCorners);

%call: STFeatureExtract
%INPUT:
%arg1: image data matrix (read in Matlab, cast to double) 
%arg2: number of channels
%arg3: [x Up Left corner, y Up Left corner, ROIwidth, ROIHeight]
%arg4: maximum number of corners to be extracted (suggested 150)
%arg5: flag enabling Harris corner detector (suggested 0)
%OUTPUT:
%arg1:column of the x coordinate of the extracted features
%arg2:column of the y coordinate of the extracted features
%arg3:number of the extracted features
%[res1,res2,res3] = STFeatureExtract(double(base), channels, [0 0 columns rows], numberCorners, flag_Harris);
%detection_points2 = [res1 res2];
%clear res1 res2 res3

%Leave only the points which are safe to be tracked:
indices = CheckPointsAndNAN(unregistered, detection_points2);
detection_points = detection_points2(indices,:);
clear indices detection_points2

%% Phase Correlation

if flag_PhaseCorrelation == 1
    base_dec = base(1:PCscaleFactor:rows,1:PCscaleFactor:columns);
    unregistered_dec = unregistered(1:PCscaleFactor:rows,1:PCscaleFactor:columns);
    [shift_xcol, shift_yrow] = ShiftByPhaseCorrelation(flag_PCglobalORlocal, base_dec(:,:,1), unregistered_dec(:,:,1));
    shift_xcol = shift_xcol*PCscaleFactor;
    shift_yrow = shift_yrow*PCscaleFactor;
    clear base_dec unregistered_dec
else
    %{
    tmp_unregistered = double(unregistered);
    tmp_unregistered(isnan(unregistered)) = nanmean(tmp_unregistered(:));
    detection_points1 = corner(tmp_unregistered, method, numberCorners);

    dx = bsxfun(@minus, detection_points(:,1), detection_points1(:,1).');
    dy = bsxfun(@minus, detection_points(:,2), detection_points1(:,2).');

    dist = dx.^2 + dy.^2;

    [val, indx1, indx2] = unique(dist(:));
    indxs = find(indx2~=[indx2(2:end); 0]);
    num = diff([indxs; length(indx2)])+1;

    keyboard
    %}

    shift_xcol = 0;
    shift_yrow = 0;
end

PointsTracked = LKTracker(double(base), double(unregistered), detection_points, [shift_xcol shift_yrow]);

%PointsTracked = [xTrackCorn yTrackCorn];
PointsBase = detection_points;

%Leave the only points inside a safe region:
indices2 = CheckPointsAndNAN(unregistered, PointsTracked);
PointsTracked = PointsTracked(indices2,:);
PointsBase = PointsBase(indices2,:);

%{
%% Features tracking

guess_xcol = detection_points(:,1)-shift_xcol;
guess_yrow = detection_points(:,2)-shift_yrow;

IndecesInsideSecondImage = find(guess_xcol>0 & guess_xcol<columns & guess_yrow>0 & guess_yrow<rows);
guess_xcol = detection_points(IndecesInsideSecondImage,1)-shift_xcol;
guess_yrow = detection_points(IndecesInsideSecondImage,2)-shift_yrow;
detection_points = detection_points(IndecesInsideSecondImage,:);
clear IndecesInsideSecondImage

%call: mexFunctionLKT
%INPUT:
%arg1: image data matrix 1 (cast to double) 
%arg2: number of channels of image 1
%arg3: image data matrix 2 (cast to double) 
%arg4: number of channels of image 2
%arg5: [xExtrCorn yExtrCorn]: 2 columns vectors --> (number of corners x 2) matrix of the extracted corners from image 1
%arg6: [xExtrCornGuess yExtrCornGuess]: 2 columns vectors --> (number of corners x 2) matrix of the "guess" corners in image 2 
%OUTPUT:
%arg1=number of features tracked reliably;
%arg2=xTrackCorn column vector of the x coords for the reliably tracked features in image 2;
%arg3=yTrackCorn column vector of the y coords for the reliably tracked features in image 2;
%arg4=xExtrTrack column vector of the x coords for the reliably extracted and tracked features in image 1;
%arg5=yExtrTrack column vector of the y coords for the reliably extracted
%and tracked features in image 1;
[numTrackFeat,xTrackCorn,yTrackCorn,xExtrTrack,yExrtTrack] = mexFunctionLKT(double(base), channels, double(unregistered), channels, detection_points, [guess_xcol guess_yrow]);
PointsTracked = [xTrackCorn yTrackCorn];
PointsBase = [xExtrTrack yExrtTrack];

%Leave the only points inside a safe region:
indices2 = CheckPointsAndNAN(unregistered, PointsTracked);
PointsTracked = PointsTracked(indices2,:);
PointsBase = PointsBase(indices2,:);
%}
