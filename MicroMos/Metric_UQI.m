function [UQI] = Metric_UQI(Frame, Reference)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% FILE NAME: Metric_UQI.m
% 
% Implementation of the UQI metric, described in the article:
% Zhou Wang and Alan C. Bovik, “A Universal Image Quality Index”, IEEE 
% Signal Processing Letters, Vol. 9, No 3, pp. 81-84, March 2002.
% The NaN values are excluded.
%
% USAGE: 
% UQI = Metric_UQI(Frame, Reference);
%
% INPUT PARAMETERS:
% 	Frame           input image of size N x M  
%   Reference       input image of size N x M 
%
% OUTPUT:
%   UQI             evaluation of the UQI metric
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
Cardinality = length(Ind);

Mean_Frame = mean(Frame(Ind));
Mean_Reference = mean(Reference(Ind));

Var_Frame = var(Frame(Ind));
Var_Reference = var(Reference(Ind));

Covariance = sum((Frame(Ind)-Mean_Frame).*(Reference(Ind)-Mean_Reference))/(Cardinality-1);
%Covariance2 = cov(Frame(Ind), Reference(Ind)); Covariance = Covariance2(1,2);
        
UQI = (4*Covariance*Mean_Frame*Mean_Reference)/((Mean_Frame.^2+Mean_Reference.^2).*(Var_Frame+Var_Reference));