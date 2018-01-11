function [Mosaic] = MosaicUpdating(Mosaic, unregistered, GLOBAL, flag_Blending, MosaicOrigin, RegistrationMode)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: MosaicUpdating
% 
% To update the current mosaic "Mosaic" stitching the input image
% "unregistered" according to the 3x3 registration matrix "GLOBAL".
%
% PARAMETERS:
%  Mosaic           Matrix of the current mosaic built stitching the images
%                   from 1 to n-1.
%  unregistered     Image n to be stitched in the current mosaic.
%  GLABAL           3x3 registration matrix to warp and stitch the image
%                   "unregisterd" into the mosaic "Mosaic".
%  flag_Blending    2 to perform a linear blending completely without
%                   seams. 1 to perform quadratic blending in the
%                   stitching regions to reduce seams in the stitching 
%                   zones. 0 = not active.
%  MosaicOrigin     x-y coordinate (x-y coordinate means the column-row
%                   coordinate) in the current "Mosaic" of the ULC (Up 
%                   Left Corner) of the first image stitched.
%  MaskOverlap      Matrix of the same size of "Mosaic". Each pixel reports
%                   the index of the last image that wrote in that specific 
%                   position.
%  InterpolationMode    Interpolation used to warp the image to be 
%                   stitched. It can be: 'bicubic' or 'bilinear'
%                   (suggested) or 'nearest'.
%  RegistrationMode     To pre-fix the registration model that must be 
%                   used to register the image "unregistered". The 
%                   registration model can be chosen between 0 (projective,
%                   suggested) or 1 (affine) or 2 (translative).
%  index            Number of images stitched on the reference image 
%                   considering also the current "unregistered". In
%                   practice, index = n-1 with n = number of images
%                   composing the mosaic.
%  Corner_Position  n-1 times [Up Left Corner, U. Right C., Low R. C., LLC,
%                   "MosaicOrigin"]. n = number of images composing the
%                   mosaic. Inside "Corner_Position" are saved the
%                   coordinated of the [ULC URC LRC LLC
%                   "MosaicOrigin"] of the single n images stitched into
%                   the mosaic. Here the coordinates are correctly scaled
%                   (+1 pixels) to find the corner coordinates inside the
%                   matrix of the mosaic. 
%
% OUTPUT:
%  Mosaic           Matrix of the mosaic built stitching the images
%                   from 1 to n.
%  MosaicOrigin     x-y coordinate (x-y coordinate means the column-row
%                   coordinate) in the current "Mosaic" of the ULC (Up 
%                   Left Corner) of the first image stitched.
%  MaskOverlap      Matrix of the same size of "Mosaic". Each pixel reports
%                   the index of the last image that wrote in that specific 
%                   position.
%  Corner_Position  n times [Up Left Corner, U. Right C., Low R. C., LLC,
%                   "MosaicOrigin"]. n = number of images composing the
%                   mosaic. Inside "Corner_Position" are saved the
%                   coordinated of the [ULC URC LRC LLC
%                   "MosaicOrigin"] of the single n images stitched into
%                   the mosaic. Here the coordinates are correctly scaled
%                   (+1 pixels) to find the corner coordinates inside the
%                   matrix of the mosaic. 
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

%% OLD MOSAIC AND NEW IMAGE WARPING IN THE SAME REFERENCE

[rowsMosaic, columnsMosaic, channelsMosaic] = size(Mosaic);
[rowsU, columnsU, channelsU] = size(unregistered);

if RegistrationMode == 2; %Translative
    modello = 'affine';
elseif RegistrationMode == 1;
    modello = 'affine';
else
    modello = 'projective';
end

unregisteredWarped = myimtransform(unregistered, modello, GLOBAL, [columnsMosaic rowsMosaic], MosaicOrigin);

clear unregistered

%% BLENDING OF THE NEW IMAGE INTO THE OLD MOSAIC

UnregisteredNotNAN  = ~isnan(unregisteredWarped); %Not NaN positions
MosaicOldNotNAN     = ~isnan(Mosaic); %Not NaN positions

no_data = ~any(MosaicOldNotNAN(:));

OverlappingPositions = (UnregisteredNotNAN & MosaicOldNotNAN);
MosaicNewPositions = (UnregisteredNotNAN & ~MosaicOldNotNAN);

clear UnregisteredNotNAN MosaicOldNotNAN

if (no_data)
  flag_Blending = 0;
end

if flag_Blending==0 %Simple stitching without blending
  Mosaic(MosaicNewPositions) = unregisteredWarped(MosaicNewPositions);

elseif flag_Blending==1 %Biquadratic blending

    Mosaic(MosaicNewPositions) = unregisteredWarped(MosaicNewPositions);
    unregisteredWarped = unregisteredWarped(OverlappingPositions);

    %Original biquadratic mask
    [X,Y] = meshgrid(0:columnsU-1,0:rowsU-1);
    BlendingOriginalMask = (1-((X-(columnsU-1)/2)./((columnsU-1)/2)).^2).*(1-((Y-(rowsU-1)/2)./((rowsU-1)/2)).^2);
    %Warped biquadratic mask
    BlendingWarpedMask = myimtransform(BlendingOriginalMask, modello, GLOBAL, [columnsMosaic rowsMosaic], MosaicOrigin);
    BlendingWarpedMask = repmat(BlendingWarpedMask, [1 1 channelsU]);
    BlendingWarpedMask = BlendingWarpedMask(OverlappingPositions);
    clear BlendingOriginalMask X Y

    Mosaic(OverlappingPositions) = (1-BlendingWarpedMask).*Mosaic(OverlappingPositions) + BlendingWarpedMask.*unregisteredWarped;

% Not optimized at all !!
elseif flag_Blending==2 || flag_Blending==3 %Blending linear with boundary values

    UnregisteredNotNANpos1  = find(isnan(unregisteredWarped(:,:,1))==0); %Not NaN positions
    MosaicOldNotNANpos1     = find(isnan(Mosaic(:,:,1))==0); %Not NaN positions

    OverlappingPositions1   = find((isnan(unregisteredWarped(:,:,1)) | (~no_data & isnan(Mosaic(:,:,1))))==0);
    MosaicNewPositions      = find(~isnan(unregisteredWarped(:,:,1)) & isnan(Mosaic(:,:,1))); 


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

            for c = 1:channelsU
                MosaicCh = Mosaic(:,:,c);
                unregisteredWarpedCh = unregisteredWarped(:,:,c);
                OverlappingRegionCh  = (1-BlendingWarpedMask(OverlappingPositionsR_indeces)).*MosaicCh(OverlappingPositionsR_indeces) + BlendingWarpedMask(OverlappingPositionsR_indeces).*unregisteredWarpedCh(OverlappingPositionsR_indeces);
                MosaicCh(OverlappingPositionsR_indeces) = OverlappingRegionCh;
                Mosaic(:,:,c) = MosaicCh;
                clear MosaicCh OverlappingRegionCh unregisteredWarpedCh MosaicCh
            end
        end
        for c = 1:channelsU
            MosaicCh = Mosaic(:,:,c);
            unregisteredWarpedCh = unregisteredWarped(:,:,c);
            MosaicCh(UnregisteredNotNANpos1_Remaining) = unregisteredWarpedCh(UnregisteredNotNANpos1_Remaining);
            Mosaic(:,:,c) = MosaicCh;
            clear MosaicCh OverlappingRegionCh unregisteredWarpedCh MosaicCh
        end
    end

end

return;
end
