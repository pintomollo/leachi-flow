function FrameOut = BleachingLUTuse(Frame, LUT)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 22 August 2013
% NAME: BleachingLUTuse
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
%  Frame            Matrix to be remapped.
%  LUT              Look-Up-Table for remapping the intensity values of 
%                   Frame. 
%
% OUTPUT:
%  FrameOut         Remapped version of Frame.
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
RemappingValMax = 255;
RemappingValMin = 0;

%% SETTINGS
Frame = double(Frame);
[row, col, ch] = size(Frame);
FrameOut = zeros(row, col);

%IndNoNaN = find(~isnan(Frame));
IndSiNaN = find(isnan(Frame));

%IntansityMin = min(LUT(:,1));
%IntansityMax = max(LUT(:,1));
NumberOfBins = length(LUT(:,1));

%% PROCESSING
% LUT remapping: values inside LUT
for i = 1:NumberOfBins
    IndI = find(Frame>LUT(i,1) & Frame<=LUT(i,2));
    FrameOut(IndI) = Frame(IndI).*LUT(i,3);
end
% LUT remapping: values below LUT
IndI = find(Frame<=LUT(1,1));
FrameOut(IndI) = Frame(IndI).*LUT(1,3);
% LUT remapping: values above LUT
IndI = find(Frame>=LUT(end,2));
FrameOut(IndI) = Frame(IndI).*LUT(end,3);
% NAN remapping
FrameOut(IndSiNaN) = NaN;

% Value above range
IndI = find(FrameOut>RemappingValMax);
FrameOut(IndI) = RemappingValMax;
% Value below range
IndI = find(FrameOut<RemappingValMin);
FrameOut(IndI) = RemappingValMin;
