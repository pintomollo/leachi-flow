function LUT = BleachingLUTbuild(Reference, Frame)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 22 August 2013
% NAME: BleachingLUTbuild
% 
% To estimate the Look-Up-Table for the photo-blaching correction.
% See details in:
% F. Piccinini, A. Bevilacqua, K. Smith and P. Horvath, Vignetting and
% photo-bleaching correction in automated fluorescence microscopy from an
% array of overlapping images. In Proceedings of the 10th IEEE
% International Symposium on Biomedical Imaging (ISBI), San Francisco, CA,
% USA, April 7-11, 2013, pp. 464 - 467.
%
% INPUT:
%  Reference        Matrix of the reference frame.
%  Frame            Matrix of the same size of Reference. The objects in 
%                   Frame have been already acquired in Reference. 
%                   Frame and Reference are perfectly overlapped images 
%                   with NaNs. 
%
% OUTPUT:
%  LUT              Look-Up-Table for remapping Frame in Reference.
%

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

%% INTERNAL PARAMETERS
NumberOfBins = 20;
MinNumberOfPointPerc = 0.02;

%% SETTINGS
if size(Frame,1)~=size(Reference,1) || size(Frame,2)~=size(Reference,2)
    error('The two input images must be of the same size')
end
Reference = double(Reference(:));
Frame = double(Frame(:));

F_noNaN = find(~isnan(Frame));
R_noNaN = find(~isnan(Reference));
Ind = intersect(F_noNaN, R_noNaN);
MinNumberOfPoint = MinNumberOfPointPerc*length(Ind);

ReferInd = Reference(Ind);
FrameInd = Frame(Ind);

IntansityMax = max(ReferInd);
IntansityMin = min(ReferInd);
IntensityBin = linspace(IntansityMin, IntansityMax, NumberOfBins);

%% PROCESSING
for i = 1:NumberOfBins-1
    LUT(i,1) = IntensityBin(i);
    LUT(i,2) = IntensityBin(i+1);
    IndI = find(ReferInd>IntensityBin(i) & ReferInd<=IntensityBin(i+1));
    if length(IndI)>MinNumberOfPoint
        ReferIndIMean = mean(ReferInd(IndI));
        FrameIndIMean = mean(FrameInd(IndI));
        LUT(i,3) = ReferIndIMean./FrameIndIMean;
    else
        LUT(i,3) = NaN;
    end
end
% NaN correction for intervals without indeces
LUT_OLD = LUT;
clear LUT
for i = 1:size(LUT_OLD,1)
    if ~isnan(LUT_OLD(i,3));
        indCopyStart = i;
        break
    end
end
for i = size(LUT_OLD,1):-1:1
    if ~isnan(LUT_OLD(i,3));
        indCopyEnd = i;
        break
    end
end
LUT = LUT_OLD(indCopyStart:indCopyEnd,:);
clear LUT_OLD
for i = 2:size(LUT,1)-1
    if isnan(LUT(i,3))
        for j = i+1:size(LUT,1)
            if ~isnan(LUT(j,3))
                LUT(i,3) = (LUT(i-1,3) + LUT(j,3))/2;
            end
        end
    end
end