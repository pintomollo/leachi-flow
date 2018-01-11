function [corners] = GriddedCorner(img, method, numberCorners, gridSize, overlap)

  if (nargin < 5)
    gridSize = [5 5];
    overlap = [0.25 0.25];
  end

  if (numel(gridSize) ~= 2)
    gridSize = gridSize([1 1]);
  end
  if (numel(overlap) ~= 2)
    overlap = overlap([1 1]);
  end

  sizes = size(img);
  blockSize = ceil(sizes ./ (gridSize - (gridSize - 1) .* overlap));
  borders = ceil((blockSize .* overlap) / 2);
  nPts = max(ceil(numberCorners*2 / prod(gridSize)), 10);
  shift = borders([2 1]) + 1;
  corners = NaN(nPts*prod(gridSize), 2);
  index = 0;

  nums = blockproc(img, blockSize, @mycorner, 'BorderSize', borders, 'UseParalle', true, 'TrimBorder', false);
  corners = unique(corners, 'rows');
  corners = corners(all(corners > 2 & bsxfun(@lt, corners, sizes([2 1])-1),2),:);

  return;

  function val = mycorner(block_struct)

    pts = corner(block_struct.data, method, nPts);
    pts = bsxfun(@plus, pts, block_struct.location([2 1]) - shift);
    val = size(pts,1);
    corners(index*nPts+[1:val],:) = pts;
    index = index + 1;

    return;
  end
end

