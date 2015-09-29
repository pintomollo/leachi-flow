% IMWARP_MEX warps an image onto a different sized matrix using the provided
% transformation matrix.
%
%   WARPED_IMG = IMWARP_MEX(IMG, MAT, IS_PROJ, WARP_SIZE, ORIGIN) warps IMG using the
%   3x3 transformation MAT onto the new WARPED_IMG of size WARP_SIZE, with its ORIGIN.
%   The warping uses a bilinear interpolation and places NaN outside of IMG.
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 05.06.2015
