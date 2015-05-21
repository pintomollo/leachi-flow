function [SNR] = Metric_SNR(Frame, Reference)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% FILE NAME: Metric_SNR.m
% 
% Implementation of the Signal to Noise Ratio (SNR) for images evaluation.
% The NaN values are excluded.
%
% USAGE: 
% SNR = Metric_SNR(Frame, Reference);
%
% INPUT PARAMETERS:
% 	Frame           input image of size N x M  
%   Reference       input image of size N x M 
%
% OUTPUT:
%   SNR             evaluation of the SNR metric
%
% See also Metric_MSE, Metric_RMSE, Metric_SNR, Metric_UQI

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

if size(Frame,1)~=size(Reference,1) || size(Frame,2)~=size(Reference,2)
    error('The two input images must be of the same size')
end

Frame = double(Frame(:));
Reference = double(Reference(:));

% TO EXCLUDE THE NaN VALUES
F_noNaN = find(~isnan(Frame));
R_noNaN = find(~isnan(Reference));
Ind = intersect(F_noNaN, R_noNaN);

SumSquareFrame = sum(Frame(Ind).^2);

SumSquareDifference = sum((Frame(Ind) - Reference(Ind)).^2);

SNR = 10*log10(SumSquareFrame/SumSquareDifference);

