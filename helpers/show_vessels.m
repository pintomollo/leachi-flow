function show_vessels(vessels)

  if (isstruct(vessels))
    vessel = [];
    if (isfield(vessels, 'border'))
      vessel = vessels;
    elseif (isfield(vessels, 'creation_params'))
      vessel = vessels.creation_params;
    end

    if (~isempty(vessel))
      figure;hold on;
      plot(vessel.center(:,1), vessel.center(:,2), 'r');
      plot(vessel.border(:,1), vessel.border(:,2), 'b');
    end
  end

  return;
end
