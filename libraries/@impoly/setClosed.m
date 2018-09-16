function setClosed(obj, is_closed)

  if (~ishandle(obj.hPolygon))
    error('impolygon graphic object deleted');
  end

  data = get(obj.hPolygon, 'userdata');

  if (is_closed ~= data.is_closed)
    if (is_closed)

      pos = data.position;
      pos = pos([1:end 1], :);

      set(obj.hEdges, 'xdata', pos(:,1), 'ydata', pos(:,2));
      set(obj.hVertex, 'xdata', pos(:,1), 'ydata', pos(:,2));
    else

      pos = data.position;
      pos = pos(1:end-1, :);

      set(obj.hEdges, 'xdata', pos(:,1), 'ydata', pos(:,2));
      set(obj.hVertex, 'xdata', pos(:,1), 'ydata', pos(:,2));
    end

    data.is_closed = is_closed;
    data.position = pos;

    set(obj.hPolygon, 'userdata', data);
  end

  return;
end
