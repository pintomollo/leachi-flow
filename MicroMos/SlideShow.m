function  SlideShow(parameters, Mosaic, MatricesGLOBAL, MosaicOrigin, LookUpTable)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: SlideShow
% 
% To save a single image for each version of the mosaic obtained during 
% the image registration stage.
%
% PARAMETERS:
%  parameters       Structure containing the parameters set into the files 
%                   "ParametersByUser.m" and "ParametersDefault.m".
%  Mosaic           Matrix of the current mosaic built stitching the images
%                   from 1 to n.
%  MatricesGLOBAL   3x3xn-1 (n = number of images composing the mosaic) 
%                   registration matrices used to warp and stitch the last 
%                   n-1 images into the mosaic.
%  MosaicOrigin     x-y coordinate (x-y coordinate means the column-row
%                   coordinate) in the current "Mosaic" of the ULC (Up 
%                   Left Corner) of the first image stitched.
%  LookUpTable      256 levels Look-Up-Table to map the grey version of the 
%                   images to be stitched into a single RGB conversion.
%                   This Look-Up-Table is used only if  
%                   "parameters.flag_LookUpTable" = 1. 
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


%% PARAMETERS SETTING

if nargin < 5
    parameters.flag_LookUpTable = 0;
end

if parameters.PixelAccuracy == 0
    fPixelAccuracy = @double;
else
    fPixelAccuracy = @single;
end

[rowsMosaic, columnsMosaic, channelsMosaic] = size(Mosaic);
%Mosaic = feval(fPixelAccuracy, Mosaic);
clear Mosaic

if parameters.RegistrationMode == 2; % Translative
    RegistrationModel = 'affine';
elseif parameters.RegistrationMode == 1; %Sarebbe l'affine
    RegistrationModel = 'affine';
else
    RegistrationModel = 'projective';
end

if parameters.flag_FlatField == 1
    DirList = dir(['VIGNETTINGFUNCTION\' '*.mat']);
    if isempty(DirList)
        parameters.flag_FlatField = 0;
    else
        Struttura = load(['VIGNETTINGFUNCTION\' DirList(1).name]);
        copyfile(['VIGNETTINGFUNCTION\' DirList(1).name], ['OUTPUT\' DirList(1).name]);
        Field = cell2mat(struct2cell(Struttura));
        Field = feval(fPixelAccuracy, Field);
        Field = Field./mean(Field(:));
        Field = Field(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);
    end
    clear Struttura DirList
end

flag_Blending = parameters.flag_Blending;

start_index = 1;
stop_index = length(parameters.ImageIndexs);

%% METRICS COMPUTATION
for index=1:stop_index

    %Image to be stitched loading and pre-processing
    strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(index));
    unregistered = imread([parameters.ImageFolder parameters.ImageBaseName strnum parameters.ImageFormat]);
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
    if parameters.flag_FlatField == 1
        [unregistered LookUpTable] = FlatFieldCorrection(unregistered, Field, parameters.flag_LookUpTable, LookUpTable);
        unregistered = feval(fPixelAccuracy, unregistered);
    end
    
    [rowsU, columnsU, channelsU] = size(unregistered);

    GLOBAL = MatricesGLOBAL(:,:,index);

    unregisteredWarped = imtransform(unregistered, maketform(RegistrationModel,GLOBAL'), parameters.InterpolationMode, ...
        'XData',[MosaicOrigin(1) columnsMosaic-1+MosaicOrigin(1)],'YData',[MosaicOrigin(2) rowsMosaic-1+MosaicOrigin(2)],...
        'XYScale',[1],...
        'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
        'fill', NaN);
    
    if index==1;
        Mosaic = unregisteredWarped;
    else
        %% BLENDING OF THE NEW IMAGE INTO THE OLD MOSAIC

        UnregisteredNotNANpos1  = find(isnan(unregisteredWarped(:,:,1))==0); %Not NaN positions
        MosaicOldNotNANpos1     = find(isnan(Mosaic(:,:,1))==0); %Not NaN positions
        OverlappingPositions1   = find(isnan(unregisteredWarped(:,:,1)+Mosaic(:,:,1))==0);
        MosaicNewPositions      = find(~isnan(unregisteredWarped(:,:,1)) & isnan(Mosaic(:,:,1))); 

        if flag_Blending==0 %Simple stitching without blending
            for c = 1:channelsMosaic
                MosaicCh = Mosaic(:,:,c);
                unregisteredWarpedCh = unregisteredWarped(:,:,c);
                MosaicCh(MosaicNewPositions) = unregisteredWarpedCh(MosaicNewPositions);
                Mosaic(:,:,c) = MosaicCh;
                clear MosaicCh unregisteredWarpedCh MosaicCh
            end
        elseif flag_Blending==1 %Biquadratic blending

            %Original biquadratic mask
            [X,Y] = meshgrid(0:columnsU-1,0:rowsU-1);
            BlendingOriginalMask = (1-((X-(columnsU-1)/2)./((columnsU-1)/2)).^2).*(1-((Y-(rowsU-1)/2)./((rowsU-1)/2)).^2);
            %Warped biquadratic mask
            BlendingWarpedMask = imtransform(double(BlendingOriginalMask), maketform(RegistrationModel,GLOBAL'), parameters.InterpolationMode, ...
                'XData',[MosaicOrigin(1) columnsMosaic-1+MosaicOrigin(1)],'YData',[MosaicOrigin(2) rowsMosaic-1+MosaicOrigin(2)],...
                'XYScale',[1],...
                'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
                'fill', NaN);
            clear BlendingOriginalMask

            for c = 1:channelsMosaic
                MosaicCh = Mosaic(:,:,c);
                unregisteredWarpedCh = unregisteredWarped(:,:,c);
                OverlappingRegionCh  = (1-BlendingWarpedMask(OverlappingPositions1)).*MosaicCh(OverlappingPositions1) + BlendingWarpedMask(OverlappingPositions1).*unregisteredWarpedCh(OverlappingPositions1);
                MosaicCh(UnregisteredNotNANpos1) = unregisteredWarpedCh(UnregisteredNotNANpos1);
                MosaicCh(OverlappingPositions1) = OverlappingRegionCh;
                Mosaic(:,:,c) = MosaicCh;
                clear MosaicCh OverlappingRegionCh unregisteredWarpedCh MosaicCh
            end

        elseif flag_Blending==2 || flag_Blending==3 %Blending linear with boundary values

            if flag_Blending==2
                BlendingComputationMode = 2;
            elseif flag_Blending==3
                BlendingComputationMode = 1;
            end

            %Original boundary
            MaskContourValue = zeros(rowsMosaic, columnsMosaic);
            MaskContourValue(OverlappingPositions1) = 1;
            boundaryOriginal = bwboundaries(MaskContourValue,'noholes');
            numberRegions = length(boundaryOriginal);
            UnregisteredNotNANpos1_Remaining = UnregisteredNotNANpos1;
            if numberRegions>=1
                for r = 1:numberRegions
                    BOriginal_yrow_xcol = boundaryOriginal{r};
                    BOriginal_indeces = sub2ind([rowsMosaic, columnsMosaic], BOriginal_yrow_xcol(:,1), BOriginal_yrow_xcol(:,2));
                    clear BOriginal_yrow_xcol
                    MaskContourValue = zeros(rowsMosaic, columnsMosaic);
                    MaskContourValue(BOriginal_indeces) = 1;
                    MaskContourValue = imfill(MaskContourValue);
                    OverlappingPositionsR_indeces = find(MaskContourValue==1);

                    %Dilated boundary
                    se = strel('disk', 1); 
                    MaskContourValue = imdilate(MaskContourValue,se);
                    boundary = bwboundaries(MaskContourValue,'noholes');
                    BOriginalDilated_yrow_xcol = boundary{1}; clear boundary
                    BOriginalDilated_indeces = sub2ind([rowsMosaic, columnsMosaic], BOriginalDilated_yrow_xcol(:,1), BOriginalDilated_yrow_xcol(:,2));
                    clear BOriginalDilated_yrow_xcol

                    %Mosaic external boundary
                    BMosaicExternal_indeces_positions = ismember(BOriginalDilated_indeces, MosaicOldNotNANpos1);
                    BMosaicExternal_indeces = BOriginalDilated_indeces(BMosaicExternal_indeces_positions);
                    clear BMosaicExternal_indeces_positions

                    %Original boundary near to mosaic external boundary
                    MaskContourValue = zeros(rowsMosaic, columnsMosaic);
                    MaskContourValue(BMosaicExternal_indeces) = 1;
                    MaskContourValue = imdilate(MaskContourValue,se);
                    BMosaicDilated_indeces = find(MaskContourValue==1);
                    BMosaicInternal_indeces_positions = ismember(BOriginal_indeces, BMosaicDilated_indeces);
                    BMosaicInternal_indeces = BOriginal_indeces(BMosaicInternal_indeces_positions);
                    clear BMosaicInternal_indeces_positions
                    clear BMosaicExternal_indeces BOriginalDilated_indeces

                    %Mask building
                    MaskContourValue = zeros(rowsMosaic, columnsMosaic)-1;
                    VOLD = 0; VNEW = 1; VROI = 2;
                    MaskContourValue(OverlappingPositionsR_indeces) = VROI;
                    MaskContourValue(BOriginal_indeces) = VNEW;
                    MaskContourValue(BMosaicInternal_indeces) = VOLD;

                    UnregisteredNotNANpos1_Remaining_positions = ismember(UnregisteredNotNANpos1_Remaining, OverlappingPositionsR_indeces);
                    UnregisteredNotNANpos1_Remaining2 = UnregisteredNotNANpos1_Remaining(UnregisteredNotNANpos1_Remaining_positions==0);
                    clear UnregisteredNotNANpos1_Remaining
                    UnregisteredNotNANpos1_Remaining = UnregisteredNotNANpos1_Remaining2;
                    clear UnregisteredNotNANpos1_Remaining2

                    BlendingWarpedMask = BlendingLinearWithContourValues(MaskContourValue, VOLD, VNEW, VROI, BlendingComputationMode);
                    clear MaskContourValue

                    for c = 1:channelsMosaic
                        MosaicCh = Mosaic(:,:,c);
                        unregisteredWarpedCh = unregisteredWarped(:,:,c);
                        OverlappingRegionCh  = (1-BlendingWarpedMask(OverlappingPositionsR_indeces)).*MosaicCh(OverlappingPositionsR_indeces) + BlendingWarpedMask(OverlappingPositionsR_indeces).*unregisteredWarpedCh(OverlappingPositionsR_indeces);
                        MosaicCh(OverlappingPositionsR_indeces) = OverlappingRegionCh;
                        Mosaic(:,:,c) = MosaicCh;
                        clear MosaicCh OverlappingRegionCh unregisteredWarpedCh MosaicCh
                    end
                end
                for c = 1:channelsMosaic
                    MosaicCh = Mosaic(:,:,c);
                    unregisteredWarpedCh = unregisteredWarped(:,:,c);
                    MosaicCh(UnregisteredNotNANpos1_Remaining) = unregisteredWarpedCh(UnregisteredNotNANpos1_Remaining);
                    Mosaic(:,:,c) = MosaicCh;
                    clear MosaicCh OverlappingRegionCh unregisteredWarpedCh MosaicCh
                end
            end

        end
        
        clear OldMosaicNaNIndexes NewFrameIndexes IndexWithoutBlending
    end

    imwrite(uint8(Mosaic), ['OUTPUT\' 'CurrentMosaic_' num2str(sprintf(parameters.NumberCharactersNumber,index)) '.png'], 'png');
    
    clear CurrentMosaic CurrentFrame
    
end