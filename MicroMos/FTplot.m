function FTplot(fig, points, type, flag_plotPointNumber)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 July 2013
% NAME: FTplot
% 
% To plot on the figbase handle the 2D features.
%
% INPUT:
%  fig              figbase handle.
%  points           points coordinates (Matrix n rows x 2 columns): 
%                   [xcol yrow].
%  type             type and colour to be used to plot the flags. E.g.:
%                   'xr' = x and red.
%  flag_plotPointNumber  1 means that near to the flag of the points are 
%                   plotted also the cardinal number of the points


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

figure(fig), hold on, 
numpoints = length(points(:,1));
for i=1:numpoints
    plot(points(i,1), points(i,2), type)
    if flag_plotPointNumber == 1
        text(points(i,1), points(i,2),num2str(i))
    end
end
hold off