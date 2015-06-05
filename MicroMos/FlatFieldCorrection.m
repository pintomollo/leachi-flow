function  [ImageO, LookUpTable] = FlatFieldCorrection(ImageI, Field, flag_LookUpTable, LookUpTable)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: FlatFieldCorrection
% 
% Each channel of the image "ImageI" is normalized pixel-wise for the 2D
% "Field" matrix. "ImageI" and "Field" must be matrices of the same x-y
% dimension.
%
% PARAMETERS:
%  unregistered     Image to be flat-field corrected.
%  Field            Vignetting function used to flat field correct the
%                   image "ImageI". "Field" and "ImageI" must be of the
%                   same x-y size.
%  flag_LookUpTable If "ImageI" is a colour image, if "flag_LookUpTable"
%                   = 1 it is modified to have a single RGB conversion for 
%                   each grey value.
%  LookUpTable      256 levels Look-Up-Table to map the grey version of the 
%                   images to be stitched into a single RGB conversion.
%                   This Look-Up-Table is used only if "flag_LookUpTable"
%                   = 1. 
%
% OUTPUT:
%  ImageO           Version flat-field corrected of the input image.
%  LookUpTable      updated version of the input Look-Up-Table.
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

[m,n,c] = size(ImageI);

if m~=size(Field,1) || n~=size(Field,2)
    error('The vignetting function must be a matrix of the same x-y dimension of the images to be corrected.');
end

if nargin < 4
    flag_LookUpTable = 0;
end

satPix = (ImageI==255);
ImageO = bsxfun(@rdivide, ImageI, Field);

if flag_LookUpTable == 1
    ImageIgrey = rgb2gray(uint8(ImageI));
    ImageOgrey = double(ImageIgrey)./Field;
    ImageOgrey(ImageIgrey==255) = 255; ImageOgrey(ImageIgrey<0) = 0;
    ImageOgrey(ImageOgrey>255)  = 255; ImageOgrey(ImageOgrey<0) = 0;

    ImageOgrey = uint8(ImageOgrey);
    [vals, indx1, indx2] = unique(ImageOgrey(:));
    vindx = vals(indx2);
    clear ImageIgrey ImageOgrey;

    news = isnan(LookUpTable(1,vals+1));
    if (any(news))
      nvals = vals(news);
      new_indxs = ismember(vindx, nvals);

      ImageO = reshape(ImageO, [m*n c]);
      [RGBs] = mymean(ImageO(new_indxs,:), 1, indx2(new_indxs));

      RGBs = [RGBs double(nvals)].';

      LookUpTable(:,nvals+1) = RGBs;
    end

    ImageO = LookUpTable(1:3,vals(indx2)+1).';
    ImageO = reshape(ImageO, [m n c]);
    ImageO(satPix) = 255;
else
    ImageO(satPix | ImageO>255) = 255;
    ImageO(ImageI <0 | ImageO <0) = 0;
end
