function BlendingMaskOut = BlendingLinearWithContourValues(MaskContourValueOrig, VOLD, VNEW, VROI, ComputationMode)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 July 2013
% NAME: BlendingLinearWithContourValues
% 
% To build a blending mask linearly interpolating the pixels of value VROI 
% between two different countour values: VOLD and VNEW. 
%
% INPUT:
%  MaskContourValueOrig Mask of the overlapping region. The boundary of the
%                   mask contains the limit values.
%  VOLD             First limit value.
%  VNEW             Second limit value, preferred in case of doubt.
%  VROI             Value of the pixels where the process must be 
%                   performed.
%  ComputationMode  the interpolation mode is define by the values 1, 2, 3
%                   1 = linear interpolation by Delaunay triangulation.
%                   2 = linear interpolation by front propagation.
%                   3 = linear interpolation using coordinates matrices.
%
% OUTPUT:
%  BlendingMaskOut  output blending mask of the same size of 
%                   MaskContourValue
%
% %Example of usage:
% VOLD = 0; %Value boundary limit low
% VNEW = 1; %Value boundary limit high
% VROI = 2; %Value ROI
% MaskContourValue = VROI*ones(100, 100); %Prima si definisce la ROI
% MaskContourValue(1,:)=VNEW; MaskContourValue(:,end)=VNEW; MaskContourValue(end,:)=VNEW; %Poi i boundary limits hight
% MaskContourValue(:,1)=VOLD; %MaskContourValue(end,:)=VOLD; %Poi i boundary limits low
% BlendingMask = BlendingLinearWithContourValues_propagation(MaskContourValue, VOLD, VNEW, VROI);

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

% Internal parameters:
if ComputationMode == 2
    PropagationEachTime = 5; %[Pixels]. The best value is 1: higher values means faster computation and less precision! 
end

%% Pre-processing MaskContourValue
[rowsOrig, columnsOrig, channelsOrig] = size(MaskContourValueOrig);
if channelsOrig > 1
    error('Blending: MaskContourValue must be a single channel matrix.')
end

[yrowVOLDOrig, xcolVOLDOrig] = find(MaskContourValueOrig==VOLD);
[yrowVNEWOrig, xcolVNEWOrig] = find(MaskContourValueOrig==VNEW);
[yrowVROIOrig, xcolVROIOrig] = find(MaskContourValueOrig==VROI);

yrowMinOrig = min([yrowVOLDOrig; yrowVNEWOrig; yrowVROIOrig]);
yrowMaxOrig = max([yrowVOLDOrig; yrowVNEWOrig; yrowVROIOrig]);
xcolMinOrig = min([xcolVOLDOrig; xcolVNEWOrig; xcolVROIOrig]);
xcolMaxOrig = max([xcolVOLDOrig; xcolVNEWOrig; xcolVROIOrig]);

%% Minimum Bounding Box
MaskContourValue = MaskContourValueOrig(yrowMinOrig:yrowMaxOrig, xcolMinOrig:xcolMaxOrig);
[rows, columns, channels] = size(MaskContourValue);
[yrowVOLD, xcolVOLD] = find(MaskContourValue==VOLD);
[yrowVNEW, xcolVNEW] = find(MaskContourValue==VNEW);
[yrowVROI, xcolVROI] = find(MaskContourValue==VROI);
indexVOLD = find(MaskContourValue==VOLD);
indexVNEW = find(MaskContourValue==VNEW);
indexVROI = find(MaskContourValue==VROI);
NumberPixelVOLD = length(yrowVOLD);
NumberPixelVNEW = length(yrowVNEW);
NumberPixelROI = length(yrowVROI);

clear MaskContourValue
VNAN = -1;
BlendingMask = VNAN*ones(rows, columns);
BlendingMask(indexVNEW) = 1;
BlendingMask(indexVOLD) = 0;

if length(indexVNEW) == 0
    BlendingMask(indexVROI) = 0;
elseif length(indexVOLD) == 0
    BlendingMask(indexVROI) = 1;
elseif length(indexVROI) ~= 0
    if ComputationMode == 1  
        [qx,qy] = meshgrid([1:1:columns], [1:1:rows]);
        qz = griddata([xcolVOLD; xcolVNEW], [yrowVOLD; yrowVNEW], [zeros(NumberPixelVOLD, 1); ones(NumberPixelVNEW, 1)], qx, qy);
        BlendingMask(indexVROI) = qz(indexVROI);
    elseif ComputationMode == 2
        %% Processing using front propagation
        disp('Blending. Please wait...')
        sizeSquare = 1 + PropagationEachTime*2;
        se = strel('square', sizeSquare); % The correct value for sizeSquare is 3! Higher values means faster computation!

        PropagationVOLD = zeros(rows, columns);
        PropagationVOLD(indexVOLD) = 1;
        CounterVOLD = zeros(rows, columns);
        CounterVOLD(indexVOLD) = 1;
        counterVOLD = 1;
        while sum(PropagationVOLD(indexVROI))/length(indexVROI)<1
            PropagationVOLD = imdilate(PropagationVOLD,se);
            updatedIndeces = find(PropagationVOLD==1 & CounterVOLD==0);
            CounterVOLD(updatedIndeces) = counterVOLD;
            counterVOLD = counterVOLD + 1;
            clear updatedIndeces
        end
        clear PropagationVOLD

        PropagationVNEW = zeros(rows, columns);
        PropagationVNEW(indexVNEW) = 1;
        CounterVNEW = zeros(rows, columns);
        CounterVNEW(indexVNEW) = 1;
        counterVNEW = 1;
        while sum(PropagationVNEW(indexVROI))/length(indexVROI)<1
            PropagationVNEW = imdilate(PropagationVNEW,se);
            updatedIndeces = find(PropagationVNEW==1 & CounterVNEW==0);
            CounterVNEW(updatedIndeces) = counterVNEW;
            counterVNEW = counterVNEW + 1;
            clear updatedIndeces
        end
        clear PropagationVNEW

        BlendingMask(indexVROI) =  CounterVOLD(indexVROI)./(CounterVOLD(indexVROI)+CounterVNEW(indexVROI));
    elseif ComputationMode == 3
        %% Processing using Matrices
        % Processing using Matrices: pixels subdivided in block to avoid memory problems.
        MaxVectorLength = 10000; %pixels %important for the memory: for images of 1000x1000 pixel use: MaxVectorLength = 10000;
        ResidualPercentageEveryTest = 0.9;

        if NumberPixelROI <= MaxVectorLength
            MaxVectorLength = NumberPixelROI-1;
        end

        flag_STOP = 1;
        while flag_STOP == 1

            BlockEnds = 1:MaxVectorLength:NumberPixelROI;
            if BlockEnds(end) ~= NumberPixelROI
                BlockEnds = [BlockEnds NumberPixelROI];
            end
            NumberBlocks = length(BlockEnds);

            try 
                for p = 1:NumberBlocks-1
                    disp(['Blending. Please wait: ' num2str(p) '/' num2str(NumberBlocks)])
                    NumberPixelROIbl = BlockEnds(p+1) - BlockEnds(p) +1;      
                    yrowVROIbl = yrowVROI(BlockEnds(p):BlockEnds(p+1));
                    xcolVROIbl = xcolVROI(BlockEnds(p):BlockEnds(p+1));
                    indexVROIbl = indexVROI(BlockEnds(p):BlockEnds(p+1));

                    MATyrowVROI_VOLD = repmat(yrowVROIbl',NumberPixelVOLD,1);
                    MATxcolVROI_VOLD = repmat(xcolVROIbl',NumberPixelVOLD,1);
                    MATyrowVROI_VNEW = repmat(yrowVROIbl',NumberPixelVNEW,1);
                    MATxcolVROI_VNEW = repmat(xcolVROIbl',NumberPixelVNEW,1);
                    clear yrowVROIbl xcolVROIbl
                    MATyrowVOLD = repmat(yrowVOLD, 1, NumberPixelROIbl);
                    MATxcolVOLD = repmat(xcolVOLD, 1, NumberPixelROIbl);
                    MATyrowVNEW = repmat(yrowVNEW, 1, NumberPixelROIbl);
                    MATxcolVNEW = repmat(xcolVNEW, 1, NumberPixelROIbl);
                    MATCurrentDistanceVOLD = sqrt((MATyrowVOLD-MATyrowVROI_VOLD).^2+(MATxcolVOLD-MATxcolVROI_VOLD).^2);
                    MATCurrentDistanceVNEW = sqrt((MATyrowVNEW-MATyrowVROI_VNEW).^2+(MATxcolVNEW-MATxcolVROI_VNEW).^2);
                    clear MATyrowVOLD MATxcolVOLD MATyrowVNEW MATxcolVNEW
                    MATMinCurrentDistanceVOLD = min(MATCurrentDistanceVOLD);
                    MATMinCurrentDistanceVNEW = min(MATCurrentDistanceVNEW);
                    MATCurrentSumDistanceVOLDVNEW = MATMinCurrentDistanceVOLD + MATMinCurrentDistanceVNEW;
                    MATCurrentValue = MATMinCurrentDistanceVOLD./MATCurrentSumDistanceVOLDVNEW;
                    clear MATMinCurrentDistanceVOLD MATMinCurrentDistanceVNEW MATCurrentSumDistanceVOLDVNEW
                    BlendingMask(indexVROIbl) = MATCurrentValue;
                end
                flag_STOP = 0;
            catch ME1
                if MaxVectorLength == 1
                    error('Blending error: MaxVectorLength=1 but processing did not succeed.')
                end
                MaxVectorLength = ceil(MaxVectorLength*ResidualPercentageEveryTest);
            end
        end
    end
end

BlendingMask(BlendingMask==VNAN) = NaN;
BlendingMaskOut = NaN*ones(rowsOrig, columnsOrig);
BlendingMaskOut(yrowMinOrig:yrowMaxOrig, xcolMinOrig:xcolMaxOrig) = BlendingMask;
clear BlendingMask
%figure, imshow(BlendingMaskOut, [], 'Border', 'Tight')

%     % Processing using Matrices: all pixels in the same time.
%     for p = 1:NumberBlocks
%         MATyrowVROI_VOLD = repmat(yrowVROI',NumberPixelVOLD,1);
%         MATxcolVROI_VOLD = repmat(xcolVROI',NumberPixelVOLD,1);
%         MATyrowVROI_VNEW = repmat(yrowVROI',NumberPixelVNEW,1);
%         MATxcolVROI_VNEW = repmat(xcolVROI',NumberPixelVNEW,1);
%         MATyrowVOLD = repmat(yrowVOLD, 1, NumberPixelROI);
%         MATxcolVOLD = repmat(xcolVOLD, 1, NumberPixelROI);
%         MATyrowVNEW = repmat(yrowVNEW, 1, NumberPixelROI);
%         MATxcolVNEW = repmat(xcolVNEW, 1, NumberPixelROI);
%         MATCurrentDistanceVOLD = sqrt((MATyrowVOLD-MATyrowVROI_VOLD).^2+(MATxcolVOLD-MATxcolVROI_VOLD).^2);
%         MATCurrentDistanceVNEW = sqrt((MATyrowVNEW-MATyrowVROI_VNEW).^2+(MATxcolVNEW-MATxcolVROI_VNEW).^2);
%         MATMinCurrentDistanceVOLD = min(MATCurrentDistanceVOLD);
%         MATMinCurrentDistanceVNEW = min(MATCurrentDistanceVNEW);
%         MATCurrentSumDistanceVOLDVNEW = MATMinCurrentDistanceVOLD + MATMinCurrentDistanceVNEW;
%         MATCurrentValue = MATMinCurrentDistanceVOLD./MATCurrentSumDistanceVOLDVNEW;
%         BlendingMask(indexVROI) = MATCurrentValue;
%     end

%     %% Processing using the "for" cycle
%     for i = 1:NumberPixelROI
%         disp(['Blending. Please wait: ' num2str(i) '/' num2str(NumberPixelROI)])
%         
%         CurrentyrowVROI = yrowVROI(i);
%         CurrentxcolVROI = xcolVROI(i);
% 
%         CurrentDistanceVOLD = sqrt((yrowVOLD-CurrentyrowVROI).^2+(xcolVOLD-CurrentxcolVROI).^2);
%         CurrentDistanceVNEW = sqrt((yrowVNEW-CurrentyrowVROI).^2+(xcolVNEW-CurrentxcolVROI).^2);
%         MinCurrentDistanceVOLD = min(CurrentDistanceVOLD);
%         MinCurrentDistanceVNEW = min(CurrentDistanceVNEW);
%         CurrentSumDistanceVOLDVNEW = MinCurrentDistanceVOLD + MinCurrentDistanceVNEW;
%         CurrentValue = MinCurrentDistanceVOLD./CurrentSumDistanceVOLDVNEW;
%         BlendingMask(CurrentyrowVROI, CurrentxcolVROI) = CurrentValue;
%     end