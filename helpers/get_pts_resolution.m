function resol = get_pts_resolution(haxes)

  if (nargin < 1)
    haxes = gca;
  elseif (~ishandle(haxes))
    resol = NaN;
    return;
  end

  switch get(haxes, 'Type')
    case 'axes'
    case 'figure'
      haxes = get(haxes, 'Child');
    otherwise
      haxes = get(haxes, 'Parent');
  end

  curr_unit = get(haxes, 'Units');
  set(haxes, 'Units', 'points');
  pos = get(haxes, 'Position');
  set(haxes, 'Units', curr_unit);

  lims = [range(get(haxes, 'XLim')) range(get(haxes, 'YLim'))];

  resol = pos(3:4) ./ lims;
  resol = mean(resol);

  return;
end
