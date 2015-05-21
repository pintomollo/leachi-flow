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

if size(ImageI,1)~=size(Field,1) || size(ImageI,2)~=size(Field,2)
    error('The vignetting function must be a matrix of the same x-y dimension of the images to be corrected.');
end

if nargin < 4
    flag_LookUpTable = 0;
end

if size(ImageI,3)==1
    %GREY correction:
    ImageO = ImageI./Field;
    ImageO(ImageI==255) = 255; ImageO(ImageI<0) = 0;
    ImageO(ImageO>255)  = 255; ImageO(ImageO<0) = 0;
else
    %RGB correction:
    % This correction strategy could generate false colours when more RGB
    % pixels of "ImageI" have the same grey conversion.
    RI = ImageI(:,:,1);
    GI = ImageI(:,:,2);
    BI = ImageI(:,:,3);

    RO = RI./Field;
    GO = GI./Field;
    BO = BI./Field;
 
    % For correcting false colours.
    if flag_LookUpTable == 1
        ImageIgrey = rgb2gray(uint8(ImageI));
        ImageOgrey = double(ImageIgrey)./double(Field);
        ImageOgrey(ImageIgrey==255) = 255; ImageOgrey(ImageIgrey<0) = 0;
        ImageOgrey(ImageOgrey>255)  = 255; ImageOgrey(ImageOgrey<0) = 0;
        ImageOgrey = uint8(ImageOgrey);
        minIOG = floor(min(ImageOgrey(:))); maxIOG = ceil(max(ImageOgrey(:)));
        for l = minIOG:maxIOG;
            positions = find(ImageOgrey==l);
            if ~isempty(positions)
                if isnan(LookUpTable(1,l+1))
                    noSatRI = find(RI(positions)~=255);
                    noSatGI = find(GI(positions)~=255);
                    noSatBI = find(BI(positions)~=255);
                    elements = intersect(intersect(noSatRI, noSatGI), noSatBI);
                    if ~isempty(elements)
                        vectorSum = RO(elements)+GO(elements)+BO(elements);
                        [vectorSumSort, vectorSumSortIndices] = sort(vectorSum);
                        value = vectorSumSortIndices(ceil(length(vectorSumSortIndices)/2));
                        position = positions(elements(value));
                    else
                        vectorSum = RO(positions)+GO(positions)+BO(positions);
                        [vectorSumSort, vectorSumSortIndices] = sort(vectorSum);
                        position = vectorSumSortIndices(ceil(length(vectorSumSortIndices)/2));
                    end
                    clear vectorSum vectorSumSort vectorSumSortIndices value elements noSatBI noSatGI noSatRI
                    LookUpTable(1,l+1) = RO(position);
                    LookUpTable(2,l+1) = GO(position);
                    LookUpTable(3,l+1) = BO(position);
                    LookUpTable(4,l+1) = l;
                end
                RO(positions) = LookUpTable(1,l+1); 
                GO(positions) = LookUpTable(2,l+1); 
                BO(positions) = LookUpTable(3,l+1); 
            end
            clear positions
        end
        clear ImageIgrey ImageOgrey minIOG maxIOG
    else
        RO(RI==255) = 255; RO(RI<0) = 0;
        GO(GI==255) = 255; GO(GI<0) = 0;
        BO(BI==255) = 255; BO(BI<0) = 0;
        RO(RO>255)  = 255; RO(RO<0) = 0;
        GO(GO>255)  = 255; GO(GO<0) = 0;
        BO(BO>255)  = 255; BO(BO<0) = 0;
    end
    
    ImageO(:,:,1) = RO; 
    ImageO(:,:,2) = GO;
    ImageO(:,:,3) = BO;
  
end