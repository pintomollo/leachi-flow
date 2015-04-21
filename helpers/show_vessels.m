function show_vessels(vessels)

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

      figure;hold on;
      plot(vessel.center(:,1), vessel.center(:,2), 'r');
      plot(vessel.border(:,1), vessel.border(:,2), 'b');
      plot(bounding_box([1 2 2 1 1]), bounding_box([3 3 4 4 3]),'k');
    end
  end

  return;
end
