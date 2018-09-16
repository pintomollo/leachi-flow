function h = impoly(varargin)
% Draws a draggable and editable polygon.

% In case we want to implement the dynamic drawing, we'll need
% Mouse motion: windowbuttonmotionfcn figure

  [reg, is_closed, contraint, color] = parseparams(varargin, 'Closed', false, 'PositionConstraintFcn', '', 'Color', 'k');

  if (isempty(reg))
    hax = gca();
    pos = [];
  elseif (numel(reg) == 1)
    if (ishghandle(reg{1}))
      hax = ancestor(reg{1}, 'axes');
      pos = [];
    else
      hax = gca ();
      pos = reg{1};
    end
  else
    if (numel(reg) > 2)
      error('Invalid number of arguments.');
    end

    hax = ancestor(reg{1}, 'axes');
    pos = reg{2};
  end

  if (isempty(hax))
    error('Not a valid axes handle.')
  end

  is_interactive = false;
  if (isempty(pos))
    pos = NaN(1,2);
    is_interactive = true;
    is_closed = false;

    hon = findall('hittest', 'on');
    hon = hon(hon ~= hax);
    set(hon, 'hittest', 'off');

  elseif (is_closed)
    if (any(pos(1,:) ~= pos(end,:)))
      pos = pos([1:end 1], :);
    end
  end

  htmp = hggroup(hax);
  addproperty('is_done', htmp, 'boolean', false);
  hlin = line(pos(:,1), pos(:,2), 'parent', htmp, 'color', color);
  hdot = line(pos(:,1), pos(:,2), 'linestyle','none','marker','+', ...
              'parent', htmp, 'color', color);

  mypoly = struct('hPolygon', htmp, ...
                  'hEdges', hlin, ...
                  'hVertex', hdot);

  mypoly = class(mypoly, 'impoly');

  polydat = struct('is_closed', is_closed, ...
                  'color', color, ...
                  'position', pos, 
                  'is_draggable', true);

  set(htmp, 'userdata', polydat);

  if (is_interactive)
    set(hdot, 'buttondownfcn', @close_polygon);

    bkp = struct('btnfnc', get(hax, 'buttondownfcn'), ...
                  'bkp', get(hax, 'userdata'), ...
                  'children', hon, ...
                  'impoly', mypoly, ...
                  'hpoly', htmp);

    set(hax, 'buttondownfcn', @create_polygon);
    set(hax, 'userdata', bkp);

    refresh(ancestor(hax, 'figure'));
    waitfor(htmp, 'is_done');
  else

    set(hlin, 'buttondownfcn', @select_edge);
    set(hdot, 'buttondownfcn', @select_vertex);
  end

  if (nargout > 0)
    h = mypoly;
  end

  return;
end

function create_polygon(hsrc, evt)

  haxe = gca;

  if (evt == 1)
    data = get(haxe, 'userdata');
    pos = get(haxe, 'currentpoint');

    prev = getPosition(data.impoly);
    if (any(isnan(prev(:))))
      prev = pos(1,1:2);
    else
      prev = [prev; pos(1,1:2)];
    end
    setPosition(data.impoly, prev);

  else

    finish_polygon(haxe)
  end

  return;
end

function close_polygon(hsrc, evt)

  haxe = gca();

  if (evt == 1)
    pos = get(haxe, 'currentpoint');

    hgrp = get(hsrc, 'parent');
    xpos = get(hsrc, 'xdata');
    ypos = get(hsrc, 'ydata');

    [dist, indx] = min((xpos - pos(1,1)) .* (xpos - pos(1,1)) + ...
                        (ypos - pos(1,2)) .* (ypos - pos(1,2)));

    if (indx == 1 && length(xpos) > 1)
      xpos = [xpos xpos(1)];
      ypos = [ypos ypos(1)];

      update_polygon(hgrp, xpos, ypos, true);
      finish_polygon(haxe);
    end    
  else

    finish_polygon(haxe);
  end

  return;
end

function select_vertex(hsrc, evt)

  haxe = gca();
  pos = get(haxe, 'currentpoint');

  hgrp = get(hsrc, 'parent');
  xpos = get(hsrc, 'xdata');
  ypos = get(hsrc, 'ydata');

  [dist, indx] = min((xpos - pos(1,1)) .* (xpos - pos(1,1)) + ...
                      (ypos - pos(1,2)) .* (ypos - pos(1,2)));

  if (evt == 1)
    hon = findobj('hittest', 'on');
    set(hon, 'hittest', 'off');
    [x, y, btn] = ginput(1);
    set(hon, 'hittest', 'on');

    if (btn == 1)
      xpos(indx) = x;
      ypos(indx) = y;

      update_polygon(hgrp, xpos, ypos);
    end
  elseif (evt == 3)

    if (indx==1)
      data = get(hgrp, 'userdata');

      if (data.is_closed)
        xpos(end) = [];
        ypos(end) = [];

        if (length(xpos)>2)
          xpos(end+1) = xpos(2);
          ypos(end+1) = ypos(2);
        end
      end
    end

    xpos(indx) = [];
    ypos(indx) = [];

    if (isempty(xpos))
      delete(hgrp);
    else
      update_polygon(hgrp, xpos, ypos);
    end
  end

  return;
end

function select_edge(hsrc, evt)

  haxe = gca();
  pos = get(haxe, 'currentpoint');

  hgrp = get(hsrc, 'parent');
  xpos = get(hsrc, 'xdata');
  ypos = get(hsrc, 'ydata');

  if (evt == 1)
    hon = findobj('hittest', 'on');
    set(hon, 'hittest', 'off');
    [x, y, btn] = ginput(1);
    set(hon, 'hittest', 'on');

    if (btn == 1)
      xpos = xpos + x - pos(1,1);
      ypos = ypos + y - pos(1,2);

      update_polygon(hgrp, xpos, ypos);
    end
  elseif (evt == 3)
    x = pos(1,1);
    y = pos(1,2);

    xvect = xpos(2:end) - xpos(1:end-1);
    yvect = ypos(2:end) - ypos(1:end-1);

    [dist, indx] = min(abs(yvect*x - xvect*y + xpos(2:end).*ypos(1:end-1) - ...
                          ypos(2:end).*xpos(1:end-1)) ./ ...
                          (yvect.*yvect + xvect.*xvect));
    xpos = [xpos(1:indx) x xpos(indx+1:end)];
    ypos = [ypos(1:indx) y ypos(indx+1:end)];

    update_polygon(hgrp, xpos, ypos);
  end

  return;
end

function finish_polygon(haxe)

  data = get(haxe, 'userdata');

  set(haxe, 'buttondownfcn', data.btnfnc);
  set(haxe, 'userdata', data.bkp);
  set(data.children,'hittest', 'on');

  childs = get(data.hpoly, 'children');
  set(childs(1), 'buttondownfcn', @select_vertex);
  set(childs(2), 'buttondownfcn', @select_edge);

  set(data.hpoly, 'is_done', true);

  return;
end

function update_polygon(hgrp, xpos, ypos, is_closed)

  set(get(hgrp, 'children'), 'xdata', xpos, 'ydata', ypos);
  data = get(hgrp, 'userdata');

  if (nargin == 4)
    data.is_closed = is_closed;
  end

  data.position = [xpos(:), ypos(:)];
  set(hgrp, 'userdata', data);

  return;
end
