function sort_mesh(nodes, links)

  figure;
  scatter(nodes(:, 1), nodes(:,2));
  hold on

  ntri = size(links, 1);

  for i=1:ntri
    triang = nodes(links(i,:),:);
    plot(triang([1:end 1],1), triang([1:end 1],2), 'k');

    keyboard;
    
    vects = triang - triang([2 3 1],:);
    lens = sqrt(sum(vects.^2, 2));
    angles = pi - acos(dot(vects, vects([2 3 1],:), 2) ./ (lens .* lens([2 3 1])));
    sperim = sum(lens)/2;
    area = sqrt(sperim*prod(sperim-lens));
    radius = prod(lens)/(2*area);
    alts = (lens .* lens([2 3 1])) / (2*radius);
    alts = alts([2 3 1]);
  end

  return;
end
