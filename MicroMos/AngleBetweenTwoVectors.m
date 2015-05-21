function [VectorAngle, VectorModule, VectorXY] = AngleBetweenTwoVectors(Vector1xy, Vector2xy)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 July 2013
% NAME: AngleBetweenTwoVectors
% 
% Given two input 2D vectors defined by their [x, y] coordinates, this
% function computes the angle, the module and the new [x, y] coordinates of
% the vector obtained performing the "vector product" operation (called
% also: cross product). 
%
% PARAMETERS:
% 	Vector1xy       Input vector coordinates (in pixels) defined as: 
%                   [Vector1x, Vector1y];
% 	Vector2xy       Input vector coordinates (in pixels) defined as:
%                   [Vector2x, Vector2y];
%
% OUTPUT:
%   VectorAngle     Degrees of direction of the output vector.
%   VectorModule    Length (in pixels) of the output vector.
% 	VectorXY        Output vector coordinates (in pixels) defined as: 
%                   [VectorX, VectorY];
%
% EXAMPLE OF USAGE: 
% [VectorAngle, VectorModule, VectorXY] = AngleBetweenTwoVectors([2,0], [0,3]);

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

Vector1x = Vector1xy(1);
Vector1y = Vector1xy(2);
Vector1abs = sqrt(Vector1x.^2 + Vector1y.^2);

Vector2x = Vector2xy(1);
Vector2y = Vector2xy(2);
Vector2abs = sqrt(Vector2x.^2 + Vector2y.^2);

VectorX = Vector1x + Vector2x;
VectorY = Vector1y + Vector2y;
VectorXY = [VectorX, VectorY];
VectorModule = sqrt(VectorX.^2 + VectorY.^2);

VectorAngleCos = acosd(VectorX/VectorModule);
if sign(VectorY)>=0
    VectorAngle = VectorAngleCos;
else
    VectorAngle = 360-VectorAngleCos;
end
