function [Vector_RMSE, Vector_Indeces] = DefineBestImageToBeRegistered(parameters, IndexBase, stop_index, ImageBase, VignettingField, RGBvignettingLookUpTable, flag_PCglobalORlocal, PCscaleFactor)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 June 2013
% NAME: DefineBestImageToBeRegistered
% 
% Run this script to define the best image to be stitched with the input
% "ImageBase" image according to the input "parameters". 
%
% INPUT:
%  parameters       Structure containing the parameters set into the files
%                   "ParametersByUser.m" and "ParametersDefault.m". Please 
%                   read carefully the two files to understand which 
%                   parameters are required.
%  IndexBase        Ordinal index of the input "ImageBase" image used as
%                   reference to define the best image to be registered.
%  stop_index       Ordinal index related to the last image that is
%                   possible to check.
%  ImageBase        Image used as reference to define the best image to be 
%                   registered.
%  VignettingField  Vignetting function used to flat field correct the
%                   images. "VignettingField" and "ImageBase" must be of 
%                   the same x-y size.
%  RGBvignettingLookUpTable      256 levels Look-Up-Table to map the grey  
%                   version of the images to be stitched into a single RGB 
%                   conversion. This Look-Up-Table is used only if 
%                   "parameters.flag_LookUpTable" == 1.
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
%  Mosaic           Matrix of the mosaic built stitching the "n" input 
%                   images according to the input parameters.
%  MaskOverlap      Matrix of the same size of "Mosaic". Each pixel reports
%                   the index of the last image (between the "n" images
%                   registered into the mosaic) that wrote in that specific
%                   position.
%  MatricesGLOBAL   3x3xn-1 ("n" = number of images composing the mosaic) 
%                   registration matrices used to warp and stitch into the 
%                   mosaic the last n-1 images.

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


%% Thresholds: "Magic Numebers" used to define the good images to be registered.
CheckRMSEmaxPercentage = 0.5; % Maximum percentage [%].
CheckANGLEmaxDeviation = 5; % Maximum deviation of the direction vector in degrees (+/-value).
CheckSHIFTmin = 2; % Minimum shift in pixel.
CheckOVERLAPminPercentage = 0.15; % Minimum percentage [%].

%% Settings
index = IndexBase;

if parameters.PixelAccuracy == 0
    fPixelAccuracy = @double;
else
    fPixelAccuracy = @single;
end

[rows, columns, channels] = size(ImageBase);
ImageBase_dec = ImageBase(1:PCscaleFactor:rows,1:PCscaleFactor:columns, :);
[row_ImageBase_dec, col_ImageBase_dec, ch_ImageBase_dec] = size(ImageBase_dec);
if (parameters.flag_Color~=0)
    if ch_ImageBase_dec == 3
        ba = rgb2gray(uint8(ImageBase_dec));
        ba = feval(fPixelAccuracy, ba);
        clear ImageBase
        ImageBase_dec = ba;
    end
end

if isempty(VignettingField)
    parameters.flag_FlatField = 0;
else
    VignettingField_dec = VignettingField(1:PCscaleFactor:rows,1:PCscaleFactor:columns, :);
end

%% Search of the good images to be registered

Vector_Angles = [];
Vector_Indeces = [];
Vector_RMSE = [];
flag_STOP = 0; % when flag_STOP == 0 means that the current tested image is not a good image to be stitched with the input ImageBase.
counter = 1; % couter related to the good image to be registered finded. 
while index < stop_index && flag_STOP == 0;
    
    %Image to be stitched loading and pre-processing
    strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(index+1));
    if strcmp(parameters.ImageFormat, '.mat')
        unregistered = load([parameters.ImageFolder parameters.ImageBaseName strnum parameters.ImageFormat]);
        unregistered = cell2mat(struct2cell(unregistered));
    else
        unregistered = imread([parameters.ImageFolder parameters.ImageBaseName strnum parameters.ImageFormat]);
    end
    unregistered = unregistered(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
    if isa(unregistered, 'uint16'); unregistered  = 255.*double(unregistered)./(2^16-1); end
    if isa(unregistered, 'uint32'); unregistered  = 255.*double(unregistered)./(2^32-1); end
    if isa(unregistered, 'uint64'); unregistered  = 255.*double(unregistered)./(2^64-1); end
    if (parameters.flag_Color==0)
        if size(unregistered, 3)~=1
            unregistered = rgb2gray(unregistered);
        end
    end
    unregistered = feval(fPixelAccuracy, unregistered);

    unregistered_dec = unregistered(1:PCscaleFactor:rows,1:PCscaleFactor:columns, :);
    if parameters.flag_FlatField == 1
        [unregistered_dec RGBvignettingLookUpTable] = FlatFieldCorrection(unregistered_dec, VignettingField_dec, parameters.flag_LookUpTable, RGBvignettingLookUpTable);
        unregistered_dec = feval(fPixelAccuracy, unregistered_dec);
    end
    if (parameters.flag_Color~=0)
        un = rgb2gray(uint8(unregistered_dec));
        un = feval(fPixelAccuracy, un);
        clear unregistered
        unregistered_dec = un;
    end

    %Registration estimation based on Phase Correlation.
    [shift_x, shift_y] = ShiftByPhaseCorrelation(flag_PCglobalORlocal, ImageBase_dec(:,:,1), unregistered_dec(:,:,1));
    
    % Check on the shift: if the unregistered_dec is shifted less than CheckSHIFTmin it is not considered at all.
    if abs(shift_x) >= CheckSHIFTmin || abs(shift_y) >= CheckSHIFTmin
        GLOBAL = [1, 0, shift_x; 0, 1, shift_y; 0, 0, 1];    
        unregisteredWarped = imtransform(double(unregistered_dec), maketform('affine',GLOBAL'), 'nearest', ...
            'XData',[0 col_ImageBase_dec-1],'YData',[0 row_ImageBase_dec-1],...
            'XYScale',[1],...
            'UData',[0 col_ImageBase_dec-1],'VData',[0 row_ImageBase_dec-1],...
            'fill', NaN);    

        Vector_RMSE(counter) = Metric_RMSE(double(ImageBase_dec), unregisteredWarped);
        Vector_Indeces(counter) = index+1;
        
        % This function is used to monitor the angle of the direction vector of the shifts.
        [Angle, canc1, canc2] = AngleBetweenTwoVectors([shift_x, 0], [0, shift_y]);
        Vector_Angles(counter) = round(Angle);
        
        if counter < 2
            % This lines are used when the first image to be registered is tested.
            ReferenceIndex = [];
            ReferenceRMSE = Vector_RMSE(1);
            DeviationMax = ReferenceRMSE*CheckRMSEmaxPercentage;
        else
            % From the second image to be registered...
            
            % Check on the RMSE: if the RMSE between unregistered_dec and 
            % ImageBase_dec is higher than ReferenceRMSE, unregistered_dec 
            % is reported as bad image to be registered.
            if Vector_RMSE(counter) <= ReferenceRMSE;
                % The reference RMSE changes if is find a negative trend
                % of the RMSE.
                ReferenceRMSE = Vector_RMSE(counter);
                ReferenceIndex = Vector_Indeces(counter);
                DeviationMax = ReferenceRMSE*CheckRMSEmaxPercentage;
            elseif Vector_RMSE(counter) >= ReferenceRMSE + DeviationMax;
                flag_STOP = 1; % The research of the good images to be registered stops here. 
            end
            
            % Check on the direction vector of the shift: if the angle of
            % the direction vector related to the shift of unregistered_dec
            % is more than CheckANGLEmaxDeviation degrees of the last one, the
            % research of the good images to be registered stops here.
            Range = mod(Vector_Angles(counter-1)-CheckANGLEmaxDeviation:Vector_Angles(counter-1)+CheckANGLEmaxDeviation,361);
            if length(find(Range==Vector_Angles(counter)))<1
                flag_STOP = 1; % The research of the good images to be registered stops here.
            end
            
            % Check on the overlap: if unregistered_dec overlaps less of
            % CheckOVERLAPminPercentage ImageBase_dec, the research of the good
            % images to be registered stops here. 
            if length(unregisteredWarped(~isnan(unregisteredWarped)))<CheckOVERLAPminPercentage*row_ImageBase_dec*col_ImageBase_dec
                flag_STOP = 1; % The research of the good images to be registered stops here.
            end
            
        end
           
        counter = counter +1;
    end
    
    index = index + 1;
end
if flag_STOP == 0
    % These lines of code are performed only if the first image check is
    % the only image that can be registered.
    Vector_RMSE = [Vector_RMSE NaN];
    Vector_Indeces = [Vector_Indeces NaN];
end
    
