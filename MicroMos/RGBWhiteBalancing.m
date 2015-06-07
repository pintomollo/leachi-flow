function IMinp = RGBWhiteBalancing(IMinp, EmptyField, WBmode)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 July 2013
% NAME: RGBWhiteBalancing
% 
% Given as input "IMinp" a RGB image, this function performs the White
% Balancing (WB) of the 3 input channels (R, G, B), using as reference the
% image "EmptyField" (image empty of objects and acquired under the same
% conditions of "IMinp"), and according to the mode: "WBmode". Practically,
% the histograms of the R, G and B channels are shifted and streched to
% have the same final mean value, but mantaining constant the ration
% between the maximum and minumum values of each single histogram. If the
% input "IMinp" is a grey level image, the input "WBmode" must be a
% positive integer value. The output "IMout" is an image of the same type
% of the image in input: if "IMinp" is RGB also "IMout" is RGB, if "IMinp"
% is a grey level image also "IMout" is a grey level image. Is important
% noting that as reference image "EmptyField" can also be used directly 
% the image "IMinp".
%
% If the input "IMinp" is a RGB image, the input "WBmode" can assume one of 
% the following values: 
% a) -4: the final mean value of the histograms is the MAX of the mean 
% values of the 3 input channels; 
% b) -3: (suggested) the final mean value of the histograms is the MEDIAN 
% of the mean values of the 3 input channels; 
% c) -2: the final mean value of the histograms is the MEAN of the mean 
% values of the 3 input channels; 
% d) -1: the final mean value of the histograms is the MIN of the mean 
% values of the 3 input channels; 
% e) N>0: the final mean value of the histograms is N (positive integer 
% value).
%
% Example of usage:
% RGBout = RGBWhiteBalancing(RGBinp); % "WBmode" = -3 by default.
% RGBout = RGBWhiteBalancing(RGBinp, RGBinp, -4);
% RGBout = RGBWhiteBalancing(RGBinp, EmptyField, -4);
% GREYout = RGBWhiteBalancing(GREYinp, GREYinp, 128);
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

%% Check on the input parameters

if nargin < 2
    EmptyField = IMinp;
end

[rowsIMinp, columnsIMinp, channelsIMinp] = size(IMinp);
[rowsEF, columnsEF, channelsEF] = size(EmptyField);

if ~isempty(EmptyField) && channelsIMinp~=channelsEF
    error('The number of channels of the input "IMinp" is different of the number of channels of the input "EmptyField".')
end

if nargin < 3
    if channels==3
        WBmode = -3; % Default for 3-channel input image
    else
        WBmode = 128;  % Default for 1-channel input image
    end
end

if channelsIMinp==3
    %% 3-channel input image
    IMinp = reshape(IMinp, [rowsIMinp*columnsIMinp channelsIMinp]);

    if (isempty(EmptyField))
      avgs = nanmean(IMinp, 1);
    else
      avgs = nanmean(reshape(EmptyField, [rowsEF*columnsEF channelsEF]), 1);
    end

    % White Balancing normalizzation factors estimation
    if WBmode == -4
        NormFact = max(avgs);
    elseif WBmode == -3
        NormFact = median(avgs);
    elseif WBmode == -2
        NormFact = mean(avgs);
    elseif WBmode == -1
        NormFact = min(avgs);
    elseif WBmode > 0
        NormFact = WBmode;
    else
        error('ERROR WHITE BALANCING: 3-channel input image. "WBmode" is not an acceptable value.')
    end
    FactorR = NormFact./avgs;

    IMinp = bsxfun(@times, IMinp, FactorR);
    IMinp = reshape(IMinp, [rowsIMinp columnsIMinp channelsIMinp]);

    %{
    % Split the reference image in the 3 channels
    Ref = EmptyField(:,:,1); Gef = EmptyField(:,:,2); Bef = EmptyField(:,:,3);
    clear EmptyField
    RefNoNaNindeces = find(isnan(Ref)==0); GefNoNaNindeces = find(isnan(Gef)==0); BefNoNaNindeces = find(isnan(Bef)==0);

    % White Balancing normalizzation factors estimation
    if WBmode == -4
        NormFact = max([mean(Ref(RefNoNaNindeces)), mean(Gef(GefNoNaNindeces)), mean(Bef(BefNoNaNindeces))]);
    elseif WBmode == -3
        NormFact = median([mean(Ref(RefNoNaNindeces)), mean(Gef(GefNoNaNindeces)), mean(Bef(BefNoNaNindeces))]);
    elseif WBmode == -2
        NormFact = mean([mean(Ref(RefNoNaNindeces)), mean(Gef(GefNoNaNindeces)), mean(Bef(BefNoNaNindeces))]);
    elseif WBmode == -1
        NormFact = min([mean(Ref(RefNoNaNindeces)), mean(Gef(GefNoNaNindeces)), mean(Bef(BefNoNaNindeces))]);
    elseif WBmode > 0
        NormFact = WBmode;
    else
        error('ERROR WHITE BALANCING: 3-channel input image. "WBmode" is not an acceptable value.')
    end
    FactorR = NormFact./mean(Ref(RefNoNaNindeces)); FactorG = NormFact./mean(Gef(GefNoNaNindeces)); FactorB = NormFact./mean(Bef(BefNoNaNindeces));
    clear Ref Gef Bef RefNoNaNindeces GefNoNaNindeces BefNoNaNindeces
    
    % Output image
    Rout = IMinp(:,:,1); Gout = IMinp(:,:,2); Bout = IMinp(:,:,3); 
    clear IMinp
    RinpNoNaNindeces = find(isnan(Rout)==0); GinpNoNaNindeces = find(isnan(Gout)==0); BinpNoNaNindeces = find(isnan(Bout)==0);
    Rout(RinpNoNaNindeces) = Rout(RinpNoNaNindeces).*FactorR; Gout(GinpNoNaNindeces) = Gout(GinpNoNaNindeces).*FactorG; Bout(BinpNoNaNindeces) = Bout(BinpNoNaNindeces).*FactorB;
    clear RinpNoNaNindeces GinpNoNaNindeces BinpNoNaNindeces
    IMout(:,:,1) = Rout; IMout(:,:,2) = Gout; IMout(:,:,3) = Bout;
    clear Rout Gout Bout
    %}
    
elseif channels==1    
    %% 1-channel input image
    
    %{
    Ref = EmptyField;
    clear EmptyField
    RefNoNaNindeces = find(isnan(Ref)==0);
    
    % White Balancing normalizzation factors estimation
    if WBmode > 0
        NormFact = WBmode;
    else
        error('ERROR WHITE BALANCING: grey input image. The normalization factor can be a positive value only.')
    end
    FactorR = NormFact./mean(Ref(RefNoNaNindeces));
    clear Ref RefNoNaNindeces

    % Output image
    Rout = IMinp; 
    clear IMinp
    RinpNoNaNindeces = find(isnan(Rout)==0);
    Rout(RinpNoNaNindeces) = Rout(RinpNoNaNindeces).*FactorR;
    IMout = Rout;
    clear Rout
   %}

    % White Balancing normalizzation factors estimation
    if WBmode > 0
        NormFact = WBmode;
    else
        error('ERROR WHITE BALANCING: grey input image. The normalization factor can be a positive value only.')
    end
    if (isempty(EmptyField))
      FactorR = NormFact./nanmean(IMinp(:));
    else
      FactorR = NormFact./nanmean(EmptyField(:));
    end
    IMinp = IMinp * FactorR;
end
