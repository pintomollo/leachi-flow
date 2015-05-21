function goodPoints = CheckPointsAndNAN(ImageI, xypoints, distance)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: CheckPointsAndNAN
% 
% To check if the x-y "xypoints" coordinates (x-y means the column-row
% coordinate inside the 2D matrix) point out pixels at a minimum distance
% "distance" pixels from a NaN pixel or a border of the 2D "ImageI" matrix.
% "goodPoints" is the list of points satisfying the above condition. 
%
% PARAMETERS:
%  ImageI           Input image used to check the corner points. 
%  xypoints         2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  distance         Number of pixel used as minimum distance between the
%                   corners and NaN values or borders.
%
% OUTPUT:
%  goodPoints       vector containing the indeces of the coordinates points 
%                   at a minimum "distance" pixels from a NaN values or a
%                   border.
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

if nargin<3
    distance = 10; %Default minimum distance in pixels
end

goodPoints = [];

for i = 1:size(xypoints,1)
    
    x=round(xypoints(i,1));
    y=round(xypoints(i,2));
    
    x1=x-distance;
    x2=x+distance;
    y1=y-distance;
    y2=y+distance;
    [SY, SX] = size(ImageI);
    if (x-distance)<1
        x1=1;
    end
    if (y-distance)<1
        y1=1;
    end
    if ((x+distance)>SX)
        x2=SX;
    end
    if ((y+distance)>SY)
        y2=SY;
    end
    if isempty(find(isnan(ImageI(y1:y2,x1:x2))==1, 1))
        if ((x>distance && x<SX-distance) && (y>distance && y<SY-distance))
            goodPoints=[goodPoints;i];
        end
    end
end