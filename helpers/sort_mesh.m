function mesh = sort_mesh(mesh)

  %figure;
  %scatter(nodes(:, 1), nodes(:,2));
  %hold on
  nodes = mesh.nodes;
  links = mesh.edges;

  ntri = size(links, 1);
  props = NaN(ntri, 2);
  params = NaN(ntri, 5);

  for i=1:ntri
    triang = nodes(links(i,:),:);

    vects = triang([2 3 1],:) - triang;
    lens = sqrt(sum(vects.^2, 2));
    %angles = acos(dot(vects, vects([2 3 1],:), 2) ./ (lens .* lens([2 3 1])));
    sperim = sum(lens)/2;
    area = sqrt(sperim*prod(sperim-lens));
    alts = 2*area ./ lens([2 3 1]);

    [height, indx] = max(alts);

    rot_indx = mod(indx, 3) + 1;
    rot_angle = acos(dot([1 0], vects(rot_indx, :)) / lens(rot_indx));
    if (vects(rot_indx, 2) < 0)
      rot_angle = 2*pi - rot_angle;
    end
    rot = [cos(-rot_angle) -sin(-rot_angle); sin(-rot_angle) cos(-rot_angle)];

    new_trig = bsxfun(@minus, triang, triang(rot_indx, :));
    new_trig = (rot * new_trig.').';
    new_trig(mod(rot_indx, 3) + 1, 2) = 0;

    base = [new_trig(indx, 1), 0];
    halfs = [base(1) lens(rot_indx)-base(1)];
    [width, side_indx] = max(halfs);

    %tmp_shift = nansum(props(:,1));

    props(i,:) = [height, width];

    if (side_indx == 1)
      params(i, :) = [height 0 abs(height./halfs) (halfs(2)<0)];
    else
      params(i, :) = [0 width abs(height./halfs) (halfs(2)<0)];
    end

    %if (any(params(i,:) < 0))

    %  plot(triang([1:end 1],1), triang([1:end 1],2), 'k');
    %  plot(new_trig([1:end 1],1), new_trig([1:end 1],2) + tmp_shift, 'k');
    %  plot([0 width width  0 0]+ base(1), [0 0 height height 0] + tmp_shift, 'b')
    %  plot([0 width] + base(1), [params(i,1) - abs(params(i,3))*(-params(i,2)) params(i,1) - abs(params(i,3))*(width-params(i,2))] + tmp_shift, 'r')
    %  plot([0 width] + base(1), [params(i,1) - abs(params(i,4))*(-params(i,2)) params(i,1) - abs(params(i,4))*(width-params(i,2))] + tmp_shift, 'r')

    %  keyboard
    %end

  end

  [vals, indxs] = sort(props(:,1), 1, 'descend');
  props = props(indxs, :);

  mesh.edges = mesh.edges(indxs, :);
  mesh.nodes = mesh.nodes(indxs, :);
  mesh.sorted = [cumsum(props(:,2)) props params(indxs, :)];
  mesh.bounding_box = [vals(1, 1), mesh.sorted(end, 1)];

  return;
end
