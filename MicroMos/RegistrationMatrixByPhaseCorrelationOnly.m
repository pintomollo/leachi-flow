function HF2F = RegistrationMatrixByPhaseCorrelationOnly(base, unregistered, PCscaleFactor, flag_PCglobalORlocal);
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: RegistrationMatrixByPhaseCorrelationOnly
% 
% To estimate the warping matrix between the image base and the image
% unregistered using the Phase Correlation Algorithm only.
%
% PARAMETERS:
%  base             Input image used to extract the corner points. We
%                   suggest to cast to double and convert to mono-channel
%                   image.
%  unregistered     Input image used to track the corner points extracted
%                   from "base". We suggest to cast to double and convert 
%                   to mono-channel image.
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
%  HF2F             3x3 registration matrix to bring the image "base" to 
%                   the image "unregistered".

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

% The corner points are extracted in the "base" image and they are tracked
% in the "unregistered" image. The functions ("STFeatureExtract.mexw32" and
% "mexFunctionLKT.mexw32") used to registered the images have been built by
% importing functions of the C++ OpenCV (see file "_licenseOpenCV.txt") and
% they require that the files with extension ".dll" are already present in
% the main folder.

if ~all(size(base)==size(unregistered))
    error('The two input matrices must be of the same dimension.');
end

[rows, columns, channels] = size(base);

%% Phase Correlation
base_dec = base(1:PCscaleFactor:rows,1:PCscaleFactor:columns);
unregistered_dec = unregistered(1:PCscaleFactor:rows,1:PCscaleFactor:columns);
[shift_xcol, shift_yrow] = ShiftByPhaseCorrelation(flag_PCglobalORlocal, base_dec(:,:,1), unregistered_dec(:,:,1));
shift_xcol = shift_xcol*PCscaleFactor;
shift_yrow = shift_yrow*PCscaleFactor; 
clear base_dec unregistered_dec

HF2F = [1, 0, -shift_xcol; 0, 1, -shift_yrow; 0, 0, 1];