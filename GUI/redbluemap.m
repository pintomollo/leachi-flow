function cMap = redbluemap(N)
% REDBLUEMAP creates a red and blue colormap as defined in [1].
%
%   REDBLUEMAP(M) returns an M-by-3 matrix containing a red and blue
%   diverging color palette. M is the number of different colors in the
%   colormap. Low values are dark blue, values in the center of the map are white,
%   and high values are dark red. If M is empty, a default value of 11 will be used.
%
% References:
%   [1] Brewer, Cynthia A., 2014. http://www.ColorBrewer.org, accessed date: 07.07.14
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2014

  % As in other colormaps, allows to color existing axes
  if (nargin < 1 || isempty(N) || N == 0)
     N = size(get(gcf,'colormap'),1);
  end

  % Get the original 11-colors map from matlab
  refMap = redbluecmap;

  % Indexes for interpolating the colors
  pos = [0:size(refMap,1)-1];
  pos = pos/pos(end);

  % Indexes for the new number of colors
  new_pos = [0:N-1];
  new_pos = new_pos/max(new_pos(end), 1);

  % Direct interpolation in each RGB coordinate
  cMap = interp1(pos(:), refMap, new_pos(:));

  return;
end
