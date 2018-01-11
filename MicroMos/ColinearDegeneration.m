function colinear = ColinearDegeneration(Set1Set2Points)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: ColinearDegeneration
% 
% To check if the input points are co-linear.
%
% PARAMETERS:
%  Set1Set2Points   x-y-z coordinates points. Set1 = Set1Set2Points(1:3,:);
%                   Set2 = Set1Set2Points(4:6,:).
%
% OUTPUT:
%  colinear         1 if the points are co-linear, 0 otherwise.
%
% I thank Peter D. Kovesi that with his Matlab functions (MATLAB and Octave
% Functions for Computer Vision and Image Processing, freely available in
% the website: http://www.csse.uwa.edu.au/~pk/research/matlabfns/) has
% contributed and inspired many functions of this tool:
% P. D. Kovesi.   
% MATLAB and Octave Functions for Computer Vision and Image Processing.
% Centre for Exploration Targeting
% School of Earth and Environment
% The University of Western Australia.
% Available from:
% <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.

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

Set1 = Set1Set2Points(1:3,:);    % ESet1Set2Pointstract Set1 and Set2 from Set1Set2Points
Set2 = Set1Set2Points(4:6,:); 
    
if size(Set1Set2Points,2) == 4
    colinear = ...
    iscolinear(Set1(:,1),Set1(:,2),Set1(:,3)) | ...
    iscolinear(Set1(:,1),Set1(:,2),Set1(:,4)) | ...
    iscolinear(Set1(:,1),Set1(:,3),Set1(:,4)) | ...
    iscolinear(Set1(:,2),Set1(:,3),Set1(:,4)) | ...
    iscolinear(Set2(:,1),Set2(:,2),Set2(:,3)) | ...
    iscolinear(Set2(:,1),Set2(:,2),Set2(:,4)) | ...
    iscolinear(Set2(:,1),Set2(:,3),Set2(:,4)) | ...
    iscolinear(Set2(:,2),Set2(:,3),Set2(:,4));
elseif size(Set1Set2Points,2) == 3
    colinear = ...
    iscolinear(Set1(:,1),Set1(:,2),Set1(:,3)) | ...
    iscolinear(Set2(:,1),Set2(:,2),Set2(:,3));
elseif size(Set1Set2Points,2) == 1
    colinear = 0;
else
    error('There must be at least 1, 3 o 4 points correspondences');
end