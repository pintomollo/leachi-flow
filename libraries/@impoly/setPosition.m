function setPosition(obj, pos, ypos)

  if (~ishandle(obj.hPolygon))
    error('impolygon graphic object deleted');
  end

  if (nargin == 3)
  	pos = [pos(:), ypos(:)];
  end

  data = get(obj.hPolygon, 'userdata');

  data.position = pos;

  set(obj.hEdges, 'xdata', pos(:,1), 'ydata', pos(:,2));
  set(obj.hVertex, 'xdata', pos(:,1), 'ydata', pos(:,2));

  data.is_closed = (size(pos, 1)>1 && all(pos(1,:) == pos(end,:)));

  set(obj.hPolygon, 'userdata', data);

  return;
end
