function hfig = show_vessels(vessels)

  if (isstruct(vessels))
    vessel = [];
    opts = [];
    if (isfield(vessels, 'border'))
      vessel = vessels;
    elseif (isfield(vessels, 'creation_params'))
      opts = vessels;
      vessel = vessels.creation_params;
    end

    if (~isempty(vessel))
      if (~isempty(opts))
        bounding_box = [1 opts.image_size(1) 1 opts.image_size(2)];
      else
        bounding_box = NaN(1, 4);
      end

      centers = vessel.center;
      x1 = centers(1:3:end, 1);
      x2 = centers(2:3:end, 1);
      y1 = centers(1:3:end, 2);
      y2 = centers(2:3:end, 2);

      vects = [x2-x1 y2-y1];
      lens = 1 ./ sum(vects.^2, 2);

      speeds = bsxfun(@times, vects, sqrt(lens));

      h = figure;hold on;
      plot(centers(:,1), centers(:,2), 'r');
      quiver(x1, y1, speeds(:,1), speeds(:,2), 'r');
      plot(vessel.border(:,1), vessel.border(:,2), 'b');
      plot(bounding_box([1 2 2 1 1]), bounding_box([3 3 4 4 3]),'k');
      plot(vessel.junction.polygon(:,1), vessel.junction.polygon(:,2), 'k')
      quiver(vessel.junction.polygon(1:5:end,1), vessel.junction.polygon(1:5:end,2), vessel.junction.vector(1,:).', vessel.junction.vector(2,:).', 'k')
    end
  end

  if (nargout > 0)
    hfig = h;
  end

  return;
end
