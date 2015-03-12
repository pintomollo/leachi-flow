function get_vessel_properties(data)

  cells = (data(:, end) == 0);
  vessel = data(~cells, :);
  cells = data(cells, :);

  vessel = vessel(:, end);
  cells = [sqrt(cells(:,1) ./ pi) cells(:,2)];

  vrange = [min(vessel)-1:max(vessel)+1].';
  vedges = [vrange(1):(vrange(end)-vrange(1))/25:vrange(end)];

  crange = [min(cells, [], 1)-1; max(cells, [], 1)+1];
  cedge1 = [crange(1, 1):(crange(2, 1)-crange(1, 1))/25:crange(2, 1)];
  cedge2 = [crange(1, 2):(crange(2, 2)-crange(1, 2))/25:crange(2, 2)];

  vcount = histc(vessel, vedges);

  options = statset('Display','final', 'MaxIter', 1000);

  vgm1 = fitgmdist(vessel, 1, 'Options', options);
  pdf1 = pdf(vgm1, vrange);
  pdf1 = max(vcount)*(pdf1 / max(pdf1));
  vgm2 = fitgmdist(vessel, 2, 'Options', options);
  pdf2 = pdf(vgm2, vrange);
  pdf2 = max(vcount)*(pdf2 / max(pdf2));
  vgm3 = fitgmdist(vessel, 3, 'Options', options);
  pdf3 = pdf(vgm3, vrange);
  pdf3 = max(vcount)*(pdf3 / max(pdf3));

  cgm1 = fitgmdist(cells, 1, 'Options', options);
  idx1 = cluster(cgm1, cells);
  count11 = histc(cells(idx1==1,1), cedge1);
  count12 = histc(cells(idx1==1,2), cedge2);
  cgm2 = fitgmdist(cells, 2, 'Options', options);
  idx2 = cluster(cgm2, cells);
  count21 = [histc(cells(idx2==1,1), cedge1), histc(cells(idx2==2,1), cedge1)];
  count22 = [histc(cells(idx2==1,2), cedge2), histc(cells(idx2==2,2), cedge2)];
  cgm3 = fitgmdist(cells, 3, 'Options', options);
  idx3 = cluster(cgm3, cells);
  count31 = [histc(cells(idx3==1,1), cedge1), histc(cells(idx3==2,1), cedge1), histc(cells(idx3==3,1), cedge1)];
  count32 = [histc(cells(idx3==1,2), cedge2), histc(cells(idx3==2,2), cedge2), histc(cells(idx3==3,2), cedge2)];

  figure;hold on
  bar(vedges, vcount, 'k');
  plot(vrange, pdf1, 'r')
  plot(vrange, pdf2, 'g')
  plot(vrange, pdf3, 'b')
  legend('counts', 'GMM 1', 'GMM 2', 'GMM 3')

  figure;
  subplot(3, 3, 1);hold on
  scatter(cells(:,1), cells(:,2), 'k')
  ezcontour(@(x,y)pdf(cgm1, [x y]), [crange(:,1).', crange(:,2).']);
  subplot(3, 3, 2);hold on
  bar(cedge1, count11, 'stacked')
  subplot(3, 3, 3);hold on
  bar(cedge2, count12, 'stacked')
  subplot(3, 3, 4);hold on
  scatter(cells(:,1), cells(:,2), 'k')
  ezcontour(@(x,y)pdf(cgm2, [x y]), [crange(:,1).', crange(:,2).']);
  subplot(3, 3, 5);hold on
  bar(cedge1, count21, 'stacked')
  subplot(3, 3, 6);hold on
  bar(cedge2, count22, 'stacked')
  subplot(3, 3, 7);hold on
  scatter(cells(:,1), cells(:,2), 'k')
  ezcontour(@(x,y)pdf(cgm3, [x y]), [crange(:,1).', crange(:,2).']);
  subplot(3, 3, 8);hold on
  bar(cedge1, count31, 'stacked')
  subplot(3, 3, 9);hold on
  bar(cedge2, count32, 'stacked')

  keyboard

  return;
end
