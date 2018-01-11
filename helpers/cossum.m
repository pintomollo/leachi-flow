function [harmo, waves] = cossum(pos, params)

  ampls = params(1,:);
  period = params(2,:);
  phases = params(3,:);

  rel_pos = bsxfun(@plus, bsxfun(@times, pos(:), 2*pi ./ period), phases);

  waves = bsxfun(@times, ampls, cos(rel_pos));
  harmo = sum(waves, 2);

  if (size(pos, 1) == 1)
    harmo = harmo.';
    waves = waves.';
  end

  return;
end
