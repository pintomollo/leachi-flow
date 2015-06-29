function simul = create_vessel(simul, opts)

  % Get the statistical properties of the vessels
  b_leachi = get_struct('botrylloides_leachi');
  vessel_props = gmdistribution(b_leachi.vessel_width.mu / opts.pixel_size, b_leachi.vessel_width.sigma / opts.pixel_size, b_leachi.vessel_width.proportions);

  % The size of the drawing area
  rim = (simul.image_size - 1) * simul.outside_ridge;
  bounding_box = [1-rim(1) simul.image_size(1)+rim(1) 1-rim(2) simul.image_size(2)+rim(2)];

  % The maximal number of trials to build the vessel
  max_trials = 20;

  % The storing structures
  vessel = get_struct('vessel');
  vessels = NaN(2, 1);
  middles = NaN(2, 0);
  junctions = NaN(2, 0);
  widths = NaN(1,0);

  % Loop over the vessels
  for i=1:simul.init_params(1)

    % And trial a number of times to create it
    for j=1:max_trials
      [new_vessel, new_middle, new_junctions, new_widths] = get_vessel(bounding_box, vessels, vessel_props, simul.init_params(2));

      % Nothing got created, something is wrong
      if (~isempty(new_vessel))
        [x,y] = polybool('&', vessels(1,:), vessels(2,:), new_vessel(1,:), new_vessel(2,:));

        if (isempty(x))
          break;
        end
      end
    end

    % We cannot create a vessel, something is wrong
    if j==max_trials
      break;
    end

    % Merge the new vessel with the previous ones
    [x,y] = polybool('|', vessels(1,:), vessels(2,:), new_vessel(1,:), new_vessel(2,:));
    vessels = [x;y];
    middles = [middles new_middle NaN(2, 1)];
    junctions = [junctions new_junctions];
    widths = [widths; new_widths];
  end

  % Polarize the vessels, i.e. define an arbitrary in/out pattern
  middles = polarize_flow(middles.');

  % Cut out the outside bits of the vessel
  [vessel.center, vessel.property] = trim_centers(middles, widths, bounding_box);

  % Clean the junctions
  vessel.border = vessels.';
  vessel.junction = refine_junctions(junctions, middles, bounding_box);

  % Create a mesh to mal the vessel
  precision = min(vessel.property(:)) * 0.75;
  [p,t] = distmesh_poly(vessels.', precision, bounding_box([1 3; 2 4]));

  % Store it
  mesh = get_struct('meshing');
  mesh.nodes = p;
  mesh.edges = t;

  % Sort the mapping
  mesh = sort_mesh(mesh);

  % Store everything
  vessel.mesh = mesh;
  simul.creation_params = vessel;

  % Display the result
  %show_vessels(simul)

  return;
end

function junct = refine_junctions(junctions, centers, bbox)

  junctions = junctions.';
  njuncs = size(junctions, 1) / 5;
  goods = all(bsxfun(@ge, junctions, bbox([1 3])) & bsxfun(@le, junctions, bbox([2 4])), 2);

  rads = NaN(1, njuncs);
  vects = NaN(2, njuncs);

  for i=njuncs:-1:1
    indxs = [1:5]+(i-1)*5;

    if (~any(goods(indxs(1:4))))
      junctions(indxs, :) = [];
      rads(i) = [];
      vects(:,i) = [];
      continue;
    end

    node = junctions(indxs(1),:);
    connecting = all(bsxfun(@eq, centers, node), 2);
    pos = find(connecting);
    direction = mod(pos, 3);
    tips = centers(pos + (3-2*direction), :);

    centered = bsxfun(@minus, [tips; junctions(indxs(2:4),:)], node);
    angles = atan2(centered(:,2), centered(:,1));

    if (direction(1)==direction(2))
      aim = angles(1:2);
      away = 3;
    elseif (direction(1)==direction(3))
      aim = angles([1 3]);
      away = 2;
    else
      aim = angles(2:3);
      away = 1;
    end

    if (aim(2) < aim(1))
      aim(2) = aim(2) + 2*pi;
    end

    pts = angles(4:end);
    pts(pts < aim(1)) = pts(pts < aim(1)) + 2*pi;
    [pts, sorting] = sort(pts, 'descend');
    sorting = sorting.' + 1;

    valids = (pts < aim(2));
    if (sum(valids) > 1)
      valids = ~valids;
    end

    pos = find(valids);
    switch pos
      case 2
        polyg = indxs([1 sorting 5]);
      case 3
        polyg = indxs([1 sorting([2 3 1]) 5]);
      otherwise
        polyg = indxs([1 sorting([3 1 2]) 5]);
    end

    junctions(indxs, :) = junctions(polyg, :);
    if (~in_triangle([0 0], centered([4 10 5 11 6 12])))

      centered(4:6, :) = centered(sorting+2,:);

      outside = [segments_intersection([0 0; centered(1,:)], centered([4:6 4],:)), ...
                 segments_intersection([0 0; centered(2,:)], centered([4:6 4],:)), ...
                 segments_intersection([0 0; centered(3,:)], centered([4:6 4],:))];

      oindx = find(~any(outside));
      output = centered(oindx,:);

      len = 1 ./ sqrt(sum(output.^2));
      new_vect = output * len;

      sides = find(sum(outside, 2) == 2);
      sides = centered(mod(sides-[1 0], 3)+4, :);

      frac = sides(:,1) * new_vect(1) + ...
             sides(:,2) * new_vect(2);

      perp = [-new_vect(2) new_vect(1)];
      perps = sides(:,1) * perp(1) + sides(:,2) * perp(2);

      repla = bsxfun(@times, perp, perps);
      new_poly = [0 0; repla(2,:); sides([2 1],:); repla(1,:); NaN(1,2)];
      new_poly = new_poly([true; frac(2)<0; true; true; frac(1)<0; true],:);

      if (size(new_poly, 1) > 5)
        new_poly = new_poly(2:end,:);
      end
      junctions(end+[1:5],:) = bsxfun(@plus, new_poly, node);

      if (direction(oindx)==2)
        new_vect = -new_vect;
      end
      vects(:,end+1) = new_vect(:);

      rads(end+1) = max(sum(bsxfun(@minus,new_poly(1,:),new_poly(2:end-1,:)).^2, 2));

      %%% Identify the two side pts using the twice-crossed segment
      %%% Compute the two 90d sideways points for the center
      %%% Check if the actual points are further back on the center (as in vascular_movement)
      %%% Create a rectangle encompassing everything, either with center or not, depending on number of nodes

      %{
      figure;hold on;
      scatter(0, 0, 'r');
      scatter(centered(4:end,1), centered(4:end,2), 'b');
      plot([zeros(3,1) centered(1:3,1)].', [zeros(3,1) centered(1:3,2)].', 'g');
      scatter(centered(sorting(pos)+2,1), centered(sorting(pos)+2,2), 'k');
      plot([0 centered(away,1)], [0 centered(away,2)], 'k');
      plot([0 centered(find(~any(outside)),1)], [0 centered(find(~any(outside)),2)], 'r');
      scatter(sides(1,1), sides(1,2), 'c');
      scatter(sides(2,1), sides(2,2), 'y');
      scatter(repla(:,1), repla(:,2), 'k');
      plot(new_poly(:,1), new_poly(:,2), 'k');
      quiver(new_poly(1,1), new_poly(1,2), new_vect(1), new_vect(2), 'b');
      %}

      junctions(indxs(1), :) = 0.5*(junctions(indxs(2), :) + junctions(indxs(4), :));
      centered_new = bsxfun(@minus, junctions(indxs(2:4),:), junctions(indxs(1), :));

      rads(i) = max(sum(centered_new.^2, 2));
      vects(:,i) = centered_new(2,:) / sqrt(sum(centered_new(2,:).^2));
    else
      rads(i) = max(sum(centered(4:end,:).^2, 2));
      vects(:,i) = centered(sorting(pos)+2,:) / sqrt(sum(centered(sorting(pos)+2,:).^2));
    end
  end

  junct = get_struct('junction');

  junct.polygon = junctions;
  junct.threshold = rads;
  junct.vector = vects;

  return;
end

function [centers, widths] = trim_centers(centers, widths, bbox)

  nsegments = size(centers, 1) / 3;
  goods = all(bsxfun(@ge, centers, bbox([1 3])) & bsxfun(@le, centers, bbox([2 4])), 2);

  corners = bbox([2 2 1 1 2; 4 3 3 4 4].');
  %figure;hold on
  %plot(centers(:,1), centers(:,2), 'r')
  %plot(corners(:,1), corners(:,2), 'k')
  %scatter(centers(goods, 1), centers(goods, 2), 'b')

  %xs = [corners(1:end-1,1) corners(2:end,1)];
  %ys = [corners(1:end-1,2) corners(2:end,2)];

  for i=nsegments:-1:1
    indxs = [1:3]+(i-1)*3;

    if (all(goods(indxs(1:2))))
      continue;
    end

    segm = centers(indxs(1:2), :);

    %scatter(segm(:,1), segm(:,2), 'g');

    %{
    origin = segm(1,:);

    segm = segm(2,:) - origin;
    len = sqrt(sum(segm.^2));

    tmp_x = xs - origin(1);
    tmp_y = ys - origin(2);

    angles = segm / len;

    rot_x = tmp_x*angles(1) + tmp_y*angles(2);
    rot_y = tmp_y*angles(1) - tmp_x*angles(2);

    dists = rot_x(:,2) + (rot_x(:,1) - rot_x(:,2)) .* rot_y(:,2) ./ (diff(rot_y, [], 2));
    valids = (dists > 0 & dists < len & xor(rot_y(:,1) > 0, rot_y(:,2) > 0));
    %}
    [valids, pos] = segments_intersection(segm, corners);

    %{
    if (any(valids))
      %pos = bsxfun(@plus, bsxfun(@times, dists(valids), angles), origin);
      rep_indx = indxs(~goods(indxs(1:2)));
      pos = pos(valids, :);
      %scatter(pos(:,1), pos(:,2), 'k');

      if (size(pos, 1) > 1)
        dists = sum(bsxfun(@minus, segm(1,:), pos).^2, 2);
        if (dists(2) < dists(1))
          rep_indx = rep_indx([2 1]);
        end
      end
      centers(rep_indx, :) = pos;
    else
      centers(indxs,:) = [];
      widths(i) = [];
    end
    %}
    if (~any(valids))
      centers(indxs,:) = [];
      widths(i) = [];
    end
  end

  %plot(centers(:,1), centers(:,2), 'b')

  %keyboard

  return;
end

function [vessel, middle, junctions, widths] = get_vessel(bounding_box, checks, props, nbranching)

  shifted = bounding_box([2 4]) + bounding_box([1 3]);
  %len = sqrt(2*sum(shifted.^2));
  len = sqrt(2*sum(shifted.^2));
  angle = rand(1)*2*pi;
  widths = random(props, 1+2*nbranching);
  center = (rand(1, 2) .* shifted - bounding_box([1 3])).';
  branches = rand(2, nbranching-1);

  corners = bounding_box([2 2 1 1 2; 4 3 3 4 4]);

  %figure;

  max_trials = 20;
  junctions = NaN(2,0);

  [middle, vessel] = get_bounding_box(len, angle, widths(1), center);

  if (nbranching == 0)
    [new_middle, new_vessel] = get_bounding_box(len, angle+pi, widths(1), center);
    middle = [middle NaN(2,1) new_middle];
    [x,y] = polybool('|', vessel(1,:), vessel(2,:), new_vessel(1,:), new_vessel(2,:));

    vessel = [x;y];

    %scatter(center(1), center(2), 'r');
  else

    orig_width = widths(1);
    orig_angle = angle;

    for i=1:nbranching

      curr_width = widths([0 1]+2*i);

      lens = [len len];
      ind = 0;
      if (i < nbranching)
        tmp_len = branches(1, i) * len/2;
        ind = 1 + (branches(2, i) > 0.5);

        lens(ind) = tmp_len;
      end

      for j=1:max_trials
        angles = randn(1, 2)*pi/2*(j/max_trials) + 2*pi/3;
        angles = orig_angle + [angles(1) -angles(2)];

        [new_middle1, new_vessel1] = get_bounding_box(lens(1), angles(1), curr_width(1), center);
        [new_middle2, new_vessel2] = get_bounding_box(lens(2), angles(2), curr_width(2), center);

        [x1,y1] = polybool('&', checks(1,:), checks(2,:), new_vessel1(1,:), new_vessel1(2,:));
        [x2,y2] = polybool('&', checks(1,:), checks(2,:), new_vessel2(1,:), new_vessel2(2,:));

        %{
        hold off;
        plot(corners(1,:), corners(2,:), 'k');
        hold on;
        plot(vessel(1,:), vessel(2,:), 'b')
        plot(checks(1,:), checks(2,:), 'c')
        plot(new_vessel1(1,:), new_vessel1(2,:), 'm')
        plot(new_vessel2(1,:), new_vessel2(2,:), 'm')
        plot([x1 x2], [y1 y2], 'r');

        drawnow
        %keyboard
        %}

        if (isempty(x1) && isempty(x2))
          break;
        end
      end

      if (j == max_trials)
        vessel = NaN(2,0);
        return;
      end

      if (ind == 1)
        [x,y] = polybool('|', vessel(1,:), vessel(2,:), new_vessel2(1,:), new_vessel2(2,:));
        checks = [x;y];
      else
        [x,y] = polybool('|', vessel(1,:), vessel(2,:), new_vessel1(1,:), new_vessel1(2,:));
        checks = [x;y];
      end

      middle = [middle NaN(2,1) new_middle1];
      [x,y] = polybool('|', vessel(1,:), vessel(2,:), new_vessel1(1,:), new_vessel1(2,:));
      vessel = [x;y];

      middle = [middle NaN(2,1) new_middle2];
      [x,y] = polybool('|', vessel(1,:), vessel(2,:), new_vessel2(1,:), new_vessel2(2,:));
      vessel = [x;y];

      %plot(vessel(1,:), vessel(2,:), 'Color', [1 0 0]*(2*i-1)/(2*nbranching))

      [vessel,new_junction] = fix_branching(vessel, center, [orig_width curr_width.'], [orig_angle angles]);
      %if (~ispolycw(better_vessel(1,:), better_vessel(2,:)))
      %  figure;hold on;
      %  plot(better_vessel(1,:), better_vessel(2,:), 'r')
      %  plot(vessel(1,:), vessel(2,:))
      %  keyboard
      %end
      %vessel = better_vessel;

      junctions = [junctions new_junction NaN(2,1)];

      %plot(vessel(1,:), vessel(2,:), 'Color', [1 0 0]*(2*i)/(2*nbranching))
      %scatter(center(1), center(2), 'r');

      if (i < nbranching)
        orig_width = curr_width(ind);
        orig_angle = angles(ind);
        center = center + [cos(orig_angle); sin(orig_angle)]*lens(ind);
        orig_angle = pi+ orig_angle;
      end
    end
  end
  [x,y] = polybool('&', vessel(1,:), vessel(2,:), corners(1,:), corners(2,:));
  vessel = [x;y];

  %{
  plot(middle(1,:), middle(2,:), 'g')
  plot(vessel(1,:), vessel(2,:), 'b')
  scatter(vessel(1, 1), vessel(2, 1), 'b')
  scatter(vessel(1, 2), vessel(2, 2), 'k')

  corners = [img_size(1)+rim(1) 1-rim(1) 1-rim(1) img_size(1)+rim(1); ...
             img_size(2)+rim(2) img_size(2)+rim(2) 1-rim(2) 1-rim(2)];
  corners = corners(:, [1:end 1]);

  [x1,y1] = poly2cw(bbox(1,:), bbox(2,:));
  [x2,y2] = poly2cw(corners(1,:), corners(2,:));
  [x,y] = polybool('&', x1, y1, x2, y2);
  vess = [x;y];

  figure;hold on;
  plot(corners(1,:), corners(2,:), 'k');
  plot([0 dvessel(1)]*50 + c(1), [0 dvessel(2)]*50 + c(2), 'g')
  scatter(c(1), c(2), 'r');
  plot(vess(1,:), vess(2,:), 'r');
  %}

  return;
end

function [vessel, junction] = fix_branching(vessel, center, widths, angles)

  %if (nargin < 5)
  %  show = false;
  %end

  ranges = [min(vessel, [], 2) max(vessel, [], 2)];
  vals = [ranges(1, 1):ranges(1, 2)];

  prec = 1024*eps;

  tans = tan(angles);
  lens = NaN(1, 3);
  bkg = NaN(1,3);

  %if (show)
  %  hfig=figure;hold on;
  %  plot(vessel(1,:), vessel(2,:), 'k')
  %  scatter(center(1), center(2), 'r');
  %end

  for i=1:length(angles)
    angle = angles(i);
    cosa = cos(angle);
    sina = sin(angle);
    if (sina < 0)
      if (cosa < 0)
        len = -widths(i) / sin(angle - pi);
      else
        len = -widths(i) / sin(2*pi - angle);
      end
    else
      if (cosa < 0)
        len = widths(i) / sin(pi - angle);
      else
        len = widths(i) / sina;
      end
    end
    lens(i) = len;
    bkg(i) = center(2) - (center(1) - len)*tans(i);

    %if (show)
    %  plot(vals, vals*tans(i) + bkg(i), 'g');
    %end
  end

  ptsx = [-(lens(1)*tans(1) + lens(2)*tans(2))/(tans(1) - tans(2)), ...
          -(lens(2)*tans(2) + lens(3)*tans(3))/(tans(2) - tans(3)), ...
          -(lens(3)*tans(3) + lens(1)*tans(1))/(tans(3) - tans(1))] + center(1);
  ptsy = ptsx .* tans + bkg;

  [in, on] = inpolygon(ptsx, ptsy, vessel(1,:), vessel(2,:));

  dist = min(bsxfun(@minus, vessel(1,:), ptsx.').^2 + bsxfun(@minus, vessel(2, :), ptsy.').^2, [], 2);
  %bads = (dist > prec & (in & ~on).');
  bads = (dist > prec & ~xor(in,on).');

  angles(angles < 0) = angles(angles < 0) + 2*pi;
  angles(angles > 2*pi) = angles(angles > 2*pi) - 2*pi;

  %if (show)
  %  scatter(ptsx, ptsy, 'b');
  %  scatter(ptsx(bads), ptsy(bads), 'm');
  %  scatter(ptsx(in), ptsy(in), 'g');
  %end

  junction = [ptsx; ptsy];

  if (any(bads))
    new_vessel = vessel;
    ptsx = ptsx(bads);
    ptsy = ptsy(bads);
    angles = angles(bads);
    bkg = bkg(bads);

    for i=1:length(angles)
      dist = sum(bsxfun(@minus, new_vessel, center).^2, 1);

      tmp_angle = atan2(new_vessel(2,:) - bkg(i), new_vessel(1,:));
      tmp_angle(tmp_angle < 0) = tmp_angle(tmp_angle < 0) + 2*pi;

      thresh = abs([tmp_angle - angles(i); pi + tmp_angle - angles(i); -pi + tmp_angle - angles(i); 2*pi + tmp_angle - angles(i); -2*pi + tmp_angle - angles(i)]);
      aligned = (min(thresh, [], 1) < prec);

      dist(~aligned) = Inf;
      [junk, index] = min(dist);

      %if (show)
      %  if (index == 1)
      %    scatter(new_vessel(1, [1 end-1]), new_vessel(2, [1 end-1]), 'y')
      %  else
      %    scatter(new_vessel(1, [index index-1]), new_vessel(2, [index index-1]), 'y')
      %  end
      %end

      if (index == 1)
        prev = new_vessel(:,end-2);
        next = new_vessel(:,2);
        others = new_vessel(:,2:end-2);
        avg = 0.5*sum(new_vessel(:, end-1:end), 2);
      elseif (index == 2)
        prev = new_vessel(:,end-1);
        next = new_vessel(:,3);
        others = new_vessel(:,3:end-1);
        avg = 0.5*sum(new_vessel(:, index-1:index), 2);
      elseif (index == size(new_vessel, 2)-1)
        prev = new_vessel(:,index-2);
        next = new_vessel(:,index+1);
        others = new_vessel(:,1:index-2);
        avg = 0.5*sum(new_vessel(:, index-1:index), 2);
      else
        prev = new_vessel(:,index-2);
        next = new_vessel(:,index+1);
        others = new_vessel(:,[index+1:end 1:index-2]);
        avg = 0.5*sum(new_vessel(:, index-1:index), 2);
      end
      intersect = segments_intersection([prev.'; ptsx(i) ptsy(i)], others.');
      if (any(intersect))
        ptsx(i) = avg(1);
        ptsy(i) = avg(2);
      else
        intersect = segments_intersection([ptsx(i) ptsy(i); next.'], others.');

        if (any(intersect))
          ptsx(i) = avg(1);
          ptsy(i) = avg(2);
        end
      end
      %if (show)
      %  plot(others(1,:), others(2,:), 'm')
      %  plot([prev(1) ptsx(i) next(1)], [prev(2) ptsy(i) next(2)], 'c')
      %  [v1, p1] = segments_intersection([prev.'; ptsx(i) ptsy(i)], others.');
      %  [v2, p2] = segments_intersection([ptsx(i) ptsy(i); next.'], others.');

      %%%%%%%%% EMERGENCY ISSUE ::: IF CROSS SOMETHING, JUST DROP THE SECOND PTS AND BASTA !!
      %  if (any(v1) || any(v2))
      %    ptsx(i) = avg(1);
      %    ptsy(i) = avg(2);
      %  end
      %end

      new_vessel(:, index) = [ptsx(i); ptsy(i)];
      if (index == 1)
        new_vessel = new_vessel(:, [1:end-2 1]);
      elseif (index == 2)
        new_vessel = new_vessel(:, [2:end-1 2]);
      else
        new_vessel(:, index-1) = [];
      end
    end

    vessel = new_vessel;

    junction(:, bads) = [ptsx; ptsy];
  end

  junction = [center junction];

  %if (show)
  %  plot(vessel(1,:), vessel(2,:), 'r')
  %  scatter(all_pts(1,:), all_pts(2,:), 'k');
  %end

  return;
end

function [middle, bbox] = get_bounding_box(len, angle, width, center)

  dvessel = [cos(angle); sin(angle)];
  dperp = [cos(angle-pi/2); sin(angle-pi/2)]*width;

  middle = [center center+len*dvessel];

  bbox = bsxfun(@plus, [len*dvessel+dperp -dvessel+dperp -dvessel-dperp len*dvessel-dperp], center);
  bbox = bbox(:, [1:end 1]);

  return;
end
