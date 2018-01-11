function [H, inliers] = RansacProjective(x1, x2, t)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 29 March 2013
% NAME: RansacProjective
% 
% To estimate the warping matrix between the points x1 and x2 (from x1 to 
% x2) using RANSAC.
%
% PARAMETERS:
%  x1               2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  x2               2-columns vector of n x-y coordinates (x-y means the
%                   column-row coordinate): 
%                   Vector_nx2 = [xp1, yp1; ...; xpn, ypn].
%  t                Maximum reprojection error (in pixels) used inside 
%                   RANSAC.
%
% OUTPUT:
%  H                3x3 registration matrix to bring the "PointsBase" to 
%                   "PointsTracked"
%  inliers          number of correspondences for which the "H" matrix is
%                   corrected according to "RANSACerror".
%
% I thank Peter D. Kovesi that with his Matlab functions (MATLAB and Octave
% Functions for Computer Vision and Image Processing, freely available in
% the website: http://www.csse.uwa.edu.au/~pk/research/matlabfns/, see
% "_licenseMIT.txt") has contributed and inspired many functions of this
% tool:
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

flag_normalization = 0;
s = 4;  % Minimum No of points needed to fit a homography.

if ~all(size(x1)==size(x2))
    error('Data sets x1 and x2 must have the same dimension');
end

[rows,npts] = size(x1);
if rows~=2 & rows~=3
    error('x1 and x2 must have 2 or 3 rows');
end

if npts < 4
    error('Must have at least 4 points to fit homography');
end

if rows == 2    % Pad data with homogeneous scale factor of 1
    x1 = [x1; ones(1,npts)];
    x2 = [x2; ones(1,npts)];        
end

if flag_normalization == 1
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
end

fittingfn = @ProjectiveModel2D;
distfn    = @ReprojectionError;
degenfn   = @ColinearDegeneration;

% x1 and x2 are 'stacked' to create a 6xN array for ransac
[H, inliers] = RANSACmodified([x1; x2], fittingfn, distfn, degenfn, s, t);

% Now do a final least squares fit on the data points considered to
% be inliers.
H = ProjectiveModel2D(x1(:,inliers), x2(:,inliers));

% Denormalise
if flag_normalization == 1
    H = T2\H*T1;
end



function Ha = ProjectiveModel2D(varargin)
    
[x1, x2] = checkargs(varargin(:));

flag_normalization = 0;

if flag_normalization == 1
    [x1c, T1] = normalise2dpts(x1);
    [x2c, T2] = normalise2dpts(x2);
    x1c = x1c(1:2,:);
    x2c = x2c(1:2,:);
    t1(1,1) = -T1(1,3)/T1(1,1);
    t1(2,1) = -T1(2,3)/T1(1,1);
    t2(1,1) = -T2(1,3)/T2(1,1);
    t2(2,1) = -T2(2,3)/T2(1,1);
else
    x1c=x1;
    x2c=x2;
end        

% A =[x1c(1,1) x1c(2,1) 1   0 0 0                 0 0 0
%     0 0 0                 x1c(1,1) x1c(2,1) 1   0 0 0
%     0 0 0                 0 0 0                 x1c(1,1) x1c(2,1) 1
%     x1c(1,2) x1c(2,2) 1   0 0 0                 0 0 0
%     0 0 0                 x1c(1,2) x1c(2,2) 1   0 0 0
%     0 0 0                 0 0 0                 x1c(1,2) x1c(2,2) 1
%     x1c(1,3) x1c(2,3) 1   0 0 0                 0 0 0
%     0 0 0                 x1c(1,3) x1c(2,3) 1   0 0 0
%     0 0 0                 0 0 0                 x1c(1,3) x1c(2,3) 1
%     x1c(1,4) x1c(2,4) 1   0 0 0                 0 0 0
%     0 0 0                 x1c(1,4) x1c(2,4) 1   0 0 0
%     0 0 0                 0 0 0                 x1c(1,4) x1c(2,4) 1];
A = [0 0 0 0 0 0 0 0 0];
for i=1:1:length(x1c)
    A = [A
        x1c(1,i) x1c(2,i) 1   0 0 0                 0 0 0
        0 0 0                 x1c(1,i) x1c(2,i) 1   0 0 0
        0 0 0                 0 0 0                 x1c(1,i) x1c(2,i) 1];
end
A = A(2:end,:);

% b = [x2c(1,1)
%     x2c(2,1)
%     x2c(1,2)
%     x2c(2,2)
%     x2c(1,3)
%     x2c(2,3)
%     x2c(1,4)
%     x2c(2,4)];
b = 0;
for i=1:1:length(x1c)
    b = [b
        x2c(1,i)
        x2c(2,i)
        1];
end
b = b(2:end,:);

%X è il vettore delle incognite 
%X = [h11
%     h12
%     h13
%     h21
%     h22
%     h23
%     h31
%     h32];

%Ricordando che: b = A * X
%Allora: X = inv(A)*b
%Che si approssima con: X = A\b
X = A\b;

Ha = [X(1) X(2) X(3)
    X(4) X(5) X(6)
    X(7) X(8) 1];

% Denormalise
if flag_normalization == 1
    Ha = T2\Ha*T1;
end
    
    
    
function [x1, x2] = checkargs(arg);
    
if length(arg) == 2
    x1 = arg{1};
    x2 = arg{2};
    if ~all(size(x1)==size(x2))
        error('x1 and x2 must have the same size');
    elseif size(x1,1) ~= 3
        error('x1 and x2 must be 3xN');
    end
elseif length(arg) == 1
    if size(arg{1},1) ~= 6
        error('Single argument x must be 6xN');
    else
        x1 = arg{1}(1:3,:);
        x2 = arg{1}(4:6,:);
    end
else
    error('Wrong number of arguments supplied');
end