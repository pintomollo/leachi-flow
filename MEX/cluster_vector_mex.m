% CLUSTER_VECTOR_MEX clusters 2D vectors at the provided resolution
%
%   [CLUST] = CLUSTER_VECTOR_MEX(X, Y, THRESH) identifies the CLUST grouping the
%   vectors [X Y] into groups closer than THRESH pixels away from each other. Note
%   that THRESH is tested on each coordinate separately. CLUST is a vector of the
%   same size as X and Y, identifying using integer values which CLUST they belong to.
%
%   [CLUST, SIZE] = CLUSTER_VECTOR_MEX(...) returns in additon the SIZE of each
%   CLUST. SIZE thus have MAX(CLUST) elements, and SIZE(N) is the number of vectors
%   in CLUST==N.
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 01.06.2015
