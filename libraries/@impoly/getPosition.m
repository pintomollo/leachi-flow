function pos = getPosition(obj)

  if (~ishandle(obj.hPolygon))
    error('impolygon graphic object deleted');
  end

  data = get(obj.hPolygon, 'userdata');
  pos = data.position;

  return;
end
