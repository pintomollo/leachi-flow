% MEDIAN_MEX constant-time median filter in C using the implementation from [1]. Note
% that the size of the memory has to be hard-coded so if you experience drastic slow-
% down, consider adapting it in MEX/median_mex.c
%
%   MED = MEDIAN_MEX(IMG, RADIUS) applies a median filter with RADIUS on IMG. The
%   kernel of the filter will thus be a 2*RADIUS+1 square.
%
%   MED = MEDIAN_MEX(IMG, RADIUS, NITER) applies the filter NITER times iteratively.
%
%   MED = MEDIAN_MEX(IMG) utilizes the default value of RADIUS=1.
%
% References:
%   [1] S. Perreault and P. Hébert, Median Filtering in Constant Time,
%       IEEE Transactions on Image Processing, (2007)
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 16.05.2014
