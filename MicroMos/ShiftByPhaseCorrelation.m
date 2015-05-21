function [xShift, yShift] = ShiftByPhaseCorrelation(flag_PCglobalORlocal, I1, I2, xGuess, yGuess)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ShiftByPhaseCorrelation
% 
% It performs the Phase Correlation between two images of the same size for
% estimating the shift of the second image (I2) respect to the first (I1). 
% Inspired to the version of the Phase Correlation algorithm 
% explained in:
% C. D. Kuglin and D. C. Hines. “The phase correlation image alignment
% method”. In: Proc. IEEE International Conference on Cybernetics and
% Society. New York, NY, U.S.A., 1975, pp. 163–165.
%
% We impose a minimum percentage of shift between the 2 input images 
% (in both x and y directions) at the 10%. Please, see the parameter 
% ShiftPercentageThreshold.
%
% USAGE: 
% [xShift, yShift] = FP_PhaseCorrelation(0, I1, I2)
% [xShift, yShift] = FP_PhaseCorrelation(1, I1, I2)
% [xShift, yShift] = FP_PhaseCorrelation(0, I1, I2, xGuess, yGuess)
% [xShift, yShift] = FP_PhaseCorrelation(1, I1, I2, xGuess, yGuess)
%
% PARAMETERS:
%   flag_PCglobalORlocal  It can assume values 0 or 1. 0 means that the 
%                   metric used to determine the best shift inside the 
%                   Phase Correlation ALgorithm is the global RMSE 
%                   performed on the entire overlapping regions. 1 means 
%                   that the used metric is the RMSE performed using only 
%                   the pixels with highest value. 
% 	I1              Input image considered as reference.
% 	I2              Input image: it is a shifted version of I1.
%   xGuess          (Optional) Initial guess relative to the x shift 
%                   (number of columns) between the image I2 and I1 with
%                   I1 as reference. 
%   yGuess          (Optional) Initial guess relative to the y shift 
%                   (number of rows) between the image I2 and I1 with
%                   I1 as reference.
%
% OUTPUT:
%   xShift          x shift (number of columns) between the image I2 and 
%                   I1 with I1 as reference.
%   yShift          y shift (number of rows) between the image I2 and 
%                   I1 with I1 as reference.
%
% See also normxcorr2

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


% INTERNAL PARAMETERS:

% flagBestCheck: the different four quadrants are not checked if it is
% zero. It enables the wrap around.
flagBestCheck = 1;

% threshold of the pcf noise: all the value inside the pcf lower of
% pcfNoiseThreshold1*(mean(pcf(:)) + pcfNoiseThreshold2*std(pcf(:))); 
% are set to zero.
pcfNoiseThreshold1 = 1;
pcfNoiseThreshold2 = 3;

% Inside the Cost Function below there are other 2 internal parameters: 
% the minimum percentage of shift checked (in both x and y directions) 
% between the two input images, the power applied to the "cost function".

% MAIN TREAD:

% The input images must be of the same size.
% Check for the size of the input images.
[height,width,ch]       = size(I1);
[height2,width2,ch2]    = size(I2);
if height~=height2 || width~=width2 || ch~=ch2
    error('The two input images must be of the same size');
end
if ch~=1
    I1 = rgb2gray(I1);
    I2 = rgb2gray(I2);
end
clear height2 width2 ch2
I1 = double(I1);
I2 = double(I2);

% Background subtraction: it improves the computation.
I1 = I1 - min(I1(:));
I2 = I2 - min(I2(:));

% FFT of each image.
F1 = fft2(double(I1));
F2 = fft2(double(I2));

% Create phase difference matrix: it's equal to the normalized Cross 
% Correlation Spectrum
pdm = exp(1i*(angle(F1)-angle(F2)));   

% Solve for phase correlation function
pcf     = ifft2(pdm);
pcf     = abs(pcf);

% Noise subtraction: it improves the next computations.
pcfNoise = pcfNoiseThreshold1*(mean(pcf(:)) + pcfNoiseThreshold2*std(pcf(:)));
pcf(pcf<pcfNoise) = 0;

if nargin == 3
    % It goes inside here if the input arguments are only 2.
    % To check for more peaks
    
    BW = imregionalmax(pcf);
    pcf(BW == 0) = 0;
    
%     pcfMaxValues = zeros(height,width); % Contains the local maximum for every region
%     BW = zeros(height,width); % Black and White image containing in white the regions with local maximums
%     BW(pcf>0) = 1;
%     L = bwlabel(BW,4);
%     NumOfL = max(L(:));
%     for l = 1:NumOfL
%         Usage = zeros(height,width);
%         Usage(L==l) = pcf(L==l);
%         [RowL, ColL] = find(Usage == max(pcf(L==l)));
%         pcfMaxValues(RowL, ColL) = pcf(RowL, ColL);
%         clear RowL ColL
%     end
%     clear pcf
%     pcf = pcfMaxValues;
%     clear pcfMaxValues Usage L BW

    % List of local maximums
    [max_val, max_ind] = sort(pcf(:));
    max_val = max_val(end:-1:1);
    
    numPicchi = 5; % Dafault 3: Maximum number of peaks checked; %@FP
    for p = 2:numPicchi
        if max_val(1) > 10*max_val(p)
            numPicchi = p-1;
            break
        end
    end
    
elseif nargin == 5
    % It goes inside here if the input arguments are only 4.
    
    % Search of the peak into a radius from the intial guess
    radius = ceil(0.05*min(height,width)); 
    
    BWregmax = imregionalmax(pcf); % Contains all the local maximums
    [RowMax, ColMax] = find(BWregmax == 1);
    NewColMax = ColMax-1-xGuess;
    NewRowMax = RowMax-1-yGuess;
    sumMax = sqrt((NewRowMax.^2) + (NewColMax.^2));

    while isempty(find(sumMax < radius))
        % If no peaks are present, it increases the radius and searches again
        radius = radius + radius;
    end
    
    % To take only the pick with the highest intensity inside a radius.
    setPositions = find(sumMax < radius);
    max_val = 0;
    yShift = yGuess;
    xShift = xGuess;
    for q = 1:length(setPositions)
        if max_val < pcf(RowMax(setPositions(q)), ColMax(setPositions(q)))
            max_val = pcf(RowMax(setPositions(q)), ColMax(setPositions(q)));
            yShift = RowMax(setPositions(q));
            xShift = ColMax(setPositions(q));
        end
    end  
   
    numPicchi = 1; % Maximum number of peaks checked
    
end

out_xShift = zeros(1,numPicchi);
out_yShift = zeros(1,numPicchi);
out_CostFuncMinValue = zeros(1,numPicchi);
for u = 1:numPicchi
    
    if nargin == 3
        % It goes inside here if the input arguments are only 2.
        [yShiftV, xShiftV] = find(abs(pcf)==max_val(u));
        yShift = yShiftV(1);
        xShift = xShiftV(1);
    end

    shift_x = xShift -1;
    shift_y = yShift -1;

    %The phase correlation uses the Fourier transform, so the initial peak 
    %is identified only in a quadrant of the image, specifically in the 
    %left upper quadrant. To see if the peak is in one of the other three 
    %quadrant, the wrap around of the coordinates is accomplished.

    tx = shift_x;
    ty = shift_y;	
    sMax.x=tx;
    sMax.y=ty;

    % 1) y x
    sMaxTemp.x=tx;
    sMaxTemp.y=ty;

    %compute the likelihood before potential wrap-around compensation 
    CostFunc = CompareOverlappingArea(flag_PCglobalORlocal, I1, I2, sMaxTemp);
    CostFuncMinValue = CostFunc;

    if flagBestCheck ~= 0

        %2) y xwa
        sMaxTemp.y = ty;
        if tx~=0
            sMaxTemp.x = -(width)+abs(tx);
            CostFunc = CompareOverlappingArea(flag_PCglobalORlocal, I1, I2, sMaxTemp);

            if CostFunc < CostFuncMinValue
                CostFuncMinValue = CostFunc;
                sMax.x = sMaxTemp.x;
                sMax.y = sMaxTemp.y;
            end
        end

        %3) ywa x
        sMaxTemp.x = tx;
        if ty~=0
            sMaxTemp.y = -(height)+abs(ty);
            CostFunc = CompareOverlappingArea(flag_PCglobalORlocal, I1, I2, sMaxTemp);	

            if CostFunc < CostFuncMinValue
                CostFuncMinValue = CostFunc;
                sMax.x = sMaxTemp.x;
                sMax.y = sMaxTemp.y;
            end
        end

        %4) ywa xwa
        if tx~=0
            sMaxTemp.x = -(width)+abs(tx);
            if ty~=0
                sMaxTemp.y = -(height)+abs(ty);
                CostFunc = CompareOverlappingArea(flag_PCglobalORlocal, I1, I2, sMaxTemp);	

                if CostFunc < CostFuncMinValue
                    CostFuncMinValue = CostFunc;
                    sMax.x = sMaxTemp.x;
                    sMax.y = sMaxTemp.y;
                end
            end

        end

    end

    xShift = sMax.x;
    yShift = sMax.y;

    out_xShift(u) = xShift;
    out_yShift(u) = yShift;
    out_CostFuncMinValue(u) = CostFuncMinValue;

end

% To take the pick with the lower CostFuncMinValue.
if isempty(out_CostFuncMinValue)
    min_CostFuncMinValue = Inf;
else
    [min_CostFuncMinValue, pos_coordinates] = min(out_CostFuncMinValue);
end
    
% OTHER TESTS:
%
%     To take the pick with the lower absolute shift. 
%     sum = sqrt((out_xShift.^2) + (out_yShift.^2));
%     [deletemi, pos_coordinates] = min(sum);
%
%     To take the pick with the lower absolute shift between a sub set of values with not high CostFuncMinValue.
%     coefB = 1.25;
%     check_pos = find(out_CostFuncMinValue < min_CostFuncMinValue*coefB);
%     if length(check_pos) > 1
%         sum = sqrt((out_xShift.^2) + (out_yShift.^2));
%         [deletemi, pos_coordinates] = min(sum);
%     end

if min_CostFuncMinValue == Inf
    xShift = NaN;
    yShift = NaN;
else
    xShift = out_xShift(pos_coordinates);
    yShift = out_yShift(pos_coordinates);
end

end


function CostFunc = CompareOverlappingArea(flag_PCglobalORlocal, I1, I2, p)
    
    % INTERNAL PARAMETERS:
    
    if ~(flag_PCglobalORlocal==0 || flag_PCglobalORlocal==1)
        error('ERROR inside Phase Correlation Algorithm: the input "flag_PCglobalORlocal" must be 0 or 1')
    end
    
    % Minimum percentage of shift between the 2 input images (in both x and y directions)
    ShiftPercentageThreshold = 0.10; % It means 10%
    
    % Percentage of pixels analyzed.
    PercentagePixelsAnalyzed = 0.01;
    
    Dx = abs(p.x);
    Dy = abs(p.y);

    %I1 = I1 / std(I1(:)); % To normalize the image but I checked and it is not an improvment
    %I2 = I2 / std(I2(:)); % To normalize the image but I checked and it is not an improvment
    
    if(p.x>=0 && p.y>=0)
        ROI_OA2 = double(I2(1:end-Dy,1:end-Dx));
        ROI_OA1 = double(I1(1+Dy:end,1+Dx:end));
    end
    
    if(p.x<0 && p.y>=0 )
		ROI_OA2 = double(I2(1:end-Dy,1+Dx:end));
        ROI_OA1 = double(I1(1+Dy:end,1:end-Dx));
    end

    if (p.x>=0 && p.y<0 )
		ROI_OA2 = double(I2(1+Dy:end,1:end-Dx));
        ROI_OA1 = double(I1(1:end-Dy,1+Dx:end));	
    end
	
    if (p.x<0 && p.y<0 )
		ROI_OA2 = double(I2(1+Dy:end,1+Dx:end));
        ROI_OA1 = double(I1(1:end-Dy,1:end-Dx));	
    end
    
%     figure(1), imshow(ROI_OA1, [], 'Border', 'Tight')
%     figure(2), imshow(ROI_OA2, [], 'Border', 'Tight') 
    
    % Often there is a problem with wrap around: ROI_OA1 and ROI_OA2 became
    % so small that they are a line of one pixel. Furthermore, sometimes
    % those lines have a CostFunc normalized smaller of the correct shift. With
    % this control ROI of a lines are excluded.
    [height,width,ch]                   = size(I1);
    [heightROI,widthROI,chROI]          = size(ROI_OA1);
    if (heightROI*widthROI)/(height*width) >= ShiftPercentageThreshold
        
        E = abs(ROI_OA2-ROI_OA1);
        if flag_PCglobalORlocal == 0
            E = E(:);
            CostFunc = sqrt(mean(E.^2));
        elseif flag_PCglobalORlocal == 1
            
            % Analysis of image's blocks
            BlockSize = 5;
            MinNumberBlockAnalyzed = 10;
            Eblock = blkproc(E, [BlockSize BlockSize], 'mean2');
            Eblock = Eblock(1:end-1, 1:end-1, :); %To delete wrap around value
            NumberBlockAnalyzed = ceil(size(Eblock,1)*size(Eblock,2)*PercentagePixelsAnalyzed);
            if NumberBlockAnalyzed < MinNumberBlockAnalyzed
                NumberBlockAnalyzed = MinNumberBlockAnalyzed;
            end
            Esort = sort(Eblock(:)); % From lowest to highest.
            Esort = Esort(end:-1:1); % From highest to lowest.
            CostFunc = sqrt(mean(Esort(1:NumberBlockAnalyzed).^2));
            
            % Analysis of image's pixels
            %Esort = sort(E(:)); % From lowest to highest.
            %Esort = Esort(end:-1:1); % From highest to lowest.
            %NumberPixelsAnalyzed = ceil(heightROI*widthROI*PercentagePixelsAnalyzed);
            %CostFunc = sqrt(mean(Esort(1:NumberPixelsAnalyzed).^2));
        end
        
    else
        CostFunc = inf;
    end
    
end