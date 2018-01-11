function [inliers, H] = ReprojectionError(H, Set1Set2Points, RANSACerror)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ReprojectionError
% 
% To estimate how many points correspondences are inlier for the "H"
% transformation matrix considering a maximum reprojection error of
% "RANSACerror" pixels.
%
% PARAMETERS:
%  H                3x3 registration matrix 
%  Set1Set2Points   x-y-z coordinates points. Set1 = Set1Set2Points(1:3,:);
%                   Set2 = Set1Set2Points(4:6,:);
%  RANSACerror      Maximum reprojection error (in pixels) used inside 
%                   RANSAC.
%
% OUTPUT:
%  inliers          number of correspondences for which the "H" matrix is
%                   corrected according to "RANSACerror".
%
% I thank Peter D. Kovesi that with his Matlab functions (MATLAB and Octave
% Functions for Computer Vision and Image Processing, freely available in
% the website: http://www.csse.uwa.edu.au/~pk/research/matlabfns/) has
% contributed and inspired many functions of this tool:
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

Set1 = Set1Set2Points(1:3,:);   % Extract Set1 and Set2 from Set1Set2Points
Set2 = Set1Set2Points(4:6,:);    

HSet1 = H*Set1;

differencesCoordinates = abs((Set2-HSet1));
distance = sqrt(differencesCoordinates(1,:).^2 + differencesCoordinates(2,:).^2);

inliers = [];
for j=1:length(distance)
    if distance(j) < RANSACerror
        inliers = [inliers, j];
    end
end
            