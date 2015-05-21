function [H, inliers]= WarpingRegistrationMode(RegistrationMode, PointsBase, PointsTracked, RANSACerror)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: WarpingRegistrationMode
% 
% To estimate the warping matrix using RANSAC. "RegistrationMode" point out
% the registration model that must be used to estimate the homography
% matrix to trasform the "PointsBase" to "PointsTracked".
%
% PARAMETERS:
%  RegistrationMode To set the registration model. RegistrationMode = 1
%                   means projective, 2 affine and 3 translative.
%  PointsBase       2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  PointsTracked    2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  RANSACerror      Maximum reprojection error (in pixels) used inside 
%                   RANSAC.
%
% OUTPUT:
%  H                3x3 registration matrix to bring the "PointsBase" to 
%                   "PointsTracked"
%  inliers          number of correspondences for which the "H" matrix is
%                   corrected according to "RANSACerror".
%
% I thank Peter D. Kovesi that with his Matlab functions (MATLAB and Octave
% Functions for Computer Vision and Image Processing, freely available in
% the website: http://www.csse.uwa.edu.au/~pk/research/matlabfns/, see
% "_licenseMIT.txt") has contributed and inspired many functions of this
% tool:
% P. D. Kovesi.   
% MATLAB and Octave Functions for Computer Vision and Image Processing.
% Centre for Exploration Targeting
% School of Earth and Environment
% The University of Western Australia.
% Available from:
% <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
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

if RegistrationMode == 2;
    [H,inliers] = RansacTranslative(PointsBase',PointsTracked',RANSACerror);
elseif RegistrationMode == 1;
    [H,inliers] = RansacAffine(PointsBase',PointsTracked',RANSACerror);
else
    [H,inliers] = RansacProjective(PointsBase',PointsTracked',RANSACerror);
end