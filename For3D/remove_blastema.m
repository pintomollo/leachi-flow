function [img, areas] = remove_blastema(img, areas, params)

  if (isempty(params) || ~isfinite(params))
    params = 0.5;
  end

  intensity = nanmean(double(img), 2);
  %stds = nanstd(double(img), [], 2);
  goods = (intensity >= max(intensity)*params);
  goods(find(goods, 1, 'last'):end) = true;
  blastema = find(~goods, 1, 'last');

  areas(1:blastema,:) = false;
  img(~areas) = 0;

  areas(1:blastema,:) = repmat(any(areas, 1), blastema, 1);
  %areas = repmat(any(areas, 1), size(areas,1), 1);

  return;
end
