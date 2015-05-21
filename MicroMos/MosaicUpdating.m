function [Mosaic, MosaicOrigin, MaskOverlap, Corner_Position] = MosaicUpdating(Mosaic, unregistered, GLOBAL, flag_Blending, MosaicOrigin, MaskOverlap, InterpolationMode, RegistrationMode, index, Corner_Position)
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

% U = Up; D = Down; L = Left; R = Right; C = Corner. 
ULC=GLOBAL*[0;0;1];
ULC=ULC./ULC(3);
DLC=GLOBAL*[0;rowsU-1;1];
DLC=DLC./DLC(3);
DRC=GLOBAL*[columnsU-1;rowsU-1;1];
DRC=DRC./DRC(3);
URC=GLOBAL*[columnsU-1;0;1];
URC=URC./URC(3);

%Maximum Bounding Box
Xmin = (min([MosaicOrigin(1),ULC(1),DLC(1)])); %ceil
Xmax = (max([columnsMosaic-1+MosaicOrigin(1),URC(1),DRC(1)])); % floor
Ymin = (min([MosaicOrigin(2),ULC(2),URC(2)])); % ceil
Ymax = (max([rowsMosaic-1+MosaicOrigin(2),DLC(2),DRC(2)])); % floor
XminI = floor(Xmin);
XmaxI = ceil(Xmax);
YminI = floor(Ymin);
YmaxI = ceil(Ymax);

unregisteredWarped = myimtransform(unregistered, modello, GLOBAL, [XminI XmaxI], [YminI YmaxI]);
%unregisteredWarped = imtransform(double(unregistered), maketform(modello,GLOBAL'), InterpolationMode, ...
%    'XData',[XminI XmaxI],'YData',[YminI YmaxI],...
%    'XYScale',[1],...
%    'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
%    'fill', NaN);    
clear unregistered

dx = XminI - MosaicOrigin(1);
dy = YminI - MosaicOrigin(2);
ox = dx;
oy = dy;

MosaicOrigin = [XminI YminI];

Mosaic = myimtransform(Mosaic, [ox XmaxI-XminI+ox], [oy YmaxI-YminI+oy]);

%Mosaic copied in the new Bounding Box. Always "nearest" interpolation.
%Mosaic = imtransform(Mosaic, maketform('affine',eye(3)), 'nearest',...
%    'XData',[ox XmaxI-XminI+ox],'YData',[oy YmaxI-YminI+oy],...
%    'XYScale',[1],...
%    'UData',[0 columnsMosaic-1],'VData',[0 rowsMosaic-1],...
%    'fill', NaN); 

[rowsMosaic2, columnsMosaic2, channelsMosaic2] = size(Mosaic);


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
    BlendingWarpedMask = myimtransform(BlendingOriginalMask, modello, GLOBAL, [XminI XmaxI], [YminI YmaxI]);
    %BlendingWarpedMask = imtransform(double(BlendingOriginalMask), maketform(modello,GLOBAL'), InterpolationMode, ...
    %    'XData',[XminI XmaxI],'YData',[YminI YmaxI],...
    %    'XYScale',[1],...
    %    'UData',[0 columnsU-1],'VData',[0 rowsU-1],...
   %     'fill', NaN);
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
    MaskContourValue = zeros(rowsMosaic2, columnsMosaic2);
    MaskContourValue(OverlappingPositions1) = 1;
    boundaryOriginal = bwboundaries(MaskContourValue,'noholes');
    numberRegions = length(boundaryOriginal);
    UnregisteredNotNANpos1_Remaining = UnregisteredNotNANpos1;
    if numberRegions>=1
        for r = 1:numberRegions
            BOriginal_yrow_xcol = boundaryOriginal{r};
            BOriginal_indeces = sub2ind([rowsMosaic2, columnsMosaic2], BOriginal_yrow_xcol(:,1), BOriginal_yrow_xcol(:,2));
            clear BOriginal_yrow_xcol
            MaskContourValue = zeros(rowsMosaic2, columnsMosaic2);
            MaskContourValue(BOriginal_indeces) = 1;
            MaskContourValue = imfill(MaskContourValue);
            OverlappingPositionsR_indeces = find(MaskContourValue==1);

            %Dilated boundary
            se = strel('disk', 1); 
            MaskContourValue = imdilate(MaskContourValue,se);
            boundary = bwboundaries(MaskContourValue,'noholes');
            BOriginalDilated_yrow_xcol = boundary{1}; clear boundary
            BOriginalDilated_indeces = sub2ind([rowsMosaic2, columnsMosaic2], BOriginalDilated_yrow_xcol(:,1), BOriginalDilated_yrow_xcol(:,2));
            clear BOriginalDilated_yrow_xcol

            %Mosaic external boundary
            BMosaicExternal_indeces_positions = ismember(BOriginalDilated_indeces, MosaicOldNotNANpos1);
            BMosaicExternal_indeces = BOriginalDilated_indeces(BMosaicExternal_indeces_positions);
            clear BMosaicExternal_indeces_positions

            %Original boundary near to mosaic external boundary
            MaskContourValue = zeros(rowsMosaic2, columnsMosaic2);
            MaskContourValue(BMosaicExternal_indeces) = 1;
            MaskContourValue = imdilate(MaskContourValue,se);
            BMosaicDilated_indeces = find(MaskContourValue==1);
            BMosaicInternal_indeces_positions = ismember(BOriginal_indeces, BMosaicDilated_indeces);
            BMosaicInternal_indeces = BOriginal_indeces(BMosaicInternal_indeces_positions);
            clear BMosaicInternal_indeces_positions
            clear BMosaicExternal_indeces BOriginalDilated_indeces

            %Mask building
            MaskContourValue = zeros(rowsMosaic2, columnsMosaic2)-1;
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


%% UPDATING OF THE MaskOverlap

if (isempty(MaskOverlap))
    MaskOverlap = uint8(zeros(rowsMosaic2, columnsMosaic2, 1));
else
    MaskOverlap = myimtransform(MaskOverlap, [ox XmaxI-XminI+ox], [oy YmaxI-YminI+oy]);
    %MaskOverlap = imtransform(MaskOverlap, maketform(modello,eye(3)), 'nearest', ... %Visto che qui è sempre uno shift intero metto subito 'nearest' per velocizzare
    %    'XData',[ox XmaxI-XminI+ox],'YData',[oy YmaxI-YminI+oy],...
    %    'XYScale',[1],...
    %    'UData',[0 columnsMosaic-1],'VData',[0 rowsMosaic-1],...
    %    'fill', NaN);
end
MaskOverlap(OverlappingPositions1) = index; 
clear OverlappingPositions1

Corner_Position = [Corner_Position, ULC, URC, DRC, DLC, [XminI YminI 1]'];
end

function [img] = myimtransform(img, modello, GLOBAL, XData, YData)

  if (nargin < 5)
    XData = modello(1):modello(2);
    YData = GLOBAL(1):GLOBAL(2);

    goodx = (XData >= 0 & XData < size(img,2));
    goody = (YData >= 0 & YData < size(img,1));

    tmp_img = NaN(length(YData), length(XData));
    tmp_img(goody, goodx) = img;
    img = tmp_img;
  else

    [X,Y] = meshgrid(XData(1):XData(2),YData(1):YData(2));
    indxs = [X(:) Y(:) zeros(numel(X), 1)] + 1;
    %Tinv = inv(GLOBAL');

    if (modello(1)=='a')
      %Tinv(1:end-1,end) = 0;
      %Tinv(end,end) = 1;

      U1 = indxs / GLOBAL';                  % Transform in homogeneous coordinates
      inv_indxs  = U1(:,1:end-1);           % Convert homogeneous coordinates to U

      %T = maketform(modello,GLOBAL');
      %[P2] = T.inverse_fcn([X(:) Y(:)]+1, T);
    else
      U1 = indxs / GLOBAL';                  % Transform in homogeneous coordinates
      inv_indxs  = bsxfun(@rdivide, U1(:,1:end-1), U1(:,end));     % Convert homogeneous coordinates to U
    end

    img = bilinear_mex(double(img), inv_indxs);

    img = reshape(img, size(X));
  end


  return;
end
