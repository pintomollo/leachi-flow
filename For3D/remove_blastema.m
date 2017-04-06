function [img, areas] = remove_blastema(img, areas, params)

  if (isempty(params) || ~isfinite(params))
    params = 0.5;
  end

  intensity = nanmean(double(img), 2);
  %[maxv,maxi,minv,mini]=local_extrema(smooth(intensity, 0.025), 40);
  %mini = [find(intensity >= max(intensity)*0.01, 1, 'first'); mini];

  %stds = nanstd(double(img), [], 2);
  goods = (smooth(intensity, 0.05) >= max(intensity)*params);
  goods(find(goods, 1, 'last'):end) = true;
  blastema = find(~goods, 1, 'last');

  %blastema2 = ceil(params*(mini(end-1)) + (1-params)*mini(end));
  %disp([blastema blastema2])
  %blastema = blastema2;

  areas(1:blastema,:) = false;
  img(~areas) = 0;

  areas(1:blastema,:) = repmat(any(areas, 1), blastema, 1);
  %areas = repmat(any(areas, 1), size(areas,1), 1);

  return;
end
