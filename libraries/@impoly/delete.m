function delete(obj)

  if (ishandle(obj.hPolygon))
    delete(obj.hPolygon);
  end

  return;
end
