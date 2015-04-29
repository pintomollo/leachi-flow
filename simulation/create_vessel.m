function opts = create_vessel(opts)

  b_leachi = get_struct('botrylloides_leachi');
  vessel_props = gmdistribution(b_leachi.vessel_width.mu / opts.pixel_size, b_leachi.vessel_width.sigma / opts.pixel_size, b_leachi.vessel_width.proportions);
  %vessel_widths = random(vessel_props, opts.init_params(1));
  %vessel_centers = randi(prod(opts.image_size), [1 opts.init_params(1)]);

  rim = (opts.image_size - 1) * opts.outside_ridge;
  bounding_box = [1-rim(1) opts.image_size(1)+rim(1) 1-rim(2) opts.image_size(2)+rim(2)];

  max_trials = 200;

  precision = min(opts.image_size)/20;

  vessel = get_struct('vessel');

  vessels = NaN(2, 1);
  middles = NaN(2, 0);
  junctions = NaN(2, 0);
  widths = NaN(1,0);
  for i=1:opts.init_params(1)
    %[centery, centerx] = ind2sub(opts.image_size, vessel_centers(i));
    %vessels(i, :) = [rand(1)*(2*pi) vessel_widths(i) centerx centery];

    %get_vessel(bounding_box, vessels(i, 1), vessels(i, 2), vessels(i,3:4), opts.init_params(2));
    for j=1:max_trials
      [new_vessel, new_middle, new_junctions, new_widths] = get_vessel(bounding_box, vessels, vessel_props, opts.init_params(2));

      if (~isempty(new_vessel))
        [x,y] = polybool('&', vessels(1,:), vessels(2,:), new_vessel(1,:), new_vessel(2,:));

        %if (i==1 && j==1)
        %figure;
        %end
        %plot(vessels(1,:), vessels(2,:), 'b')
        %hold on
        %plot(new_vessel(1,:), new_vessel(2,:), 'k')
        %plot(x, y, 'r')
        %drawnow
        %hold off

        if (isempty(x))
          break;
        end
      end
    end

    if j==max_trials
      break;
    end

    [x,y] = polybool('|', vessels(1,:), vessels(2,:), new_vessel(1,:), new_vessel(2,:));
    vessels = [x;y];
    middles = [middles new_middle NaN(2, 1)];
    junctions = [junctions new_junctions];
    widths = [widths; new_widths];

    %p1 = get_bounding_box(vessels(i,:), opts);
    %p2 = get_bounding_box(vessels(i,:)+[pi 0 0 0], opts);
    %[x,y] = polybool('|', p1(1,:), p1(2,:), p2(1,:), p2(2,:));
    %plot(x, y, 'b');
  end

  %figure;
  middles = polarize_flow(middles.');

  [vessel.center, vessel.property] = trim_centers(middles, widths, bounding_box);
  %vessel.center = middles;
  %vessel.property = widths;

  vessel.border = vessels.';
  vessel.junction = refine_junctions(junctions, middles, bounding_box);

  precision = min(vessel.property(:)) * 0.75;

  [p,t] = distmesh_poly(vessels.', precision, bounding_box([1 3; 2 4]));

  mesh = get_struct('meshing');
  mesh.nodes = p;
  mesh.edges = t;

  mesh = sort_mesh(mesh);

  vessel.mesh = mesh;

  opts.creation_params = vessel;

  %opts.creation_params = vessels(1, 1:2);
  %opts.movement_params = [cos(vessels(1,1)) sin(vessels(1, 1))];

  figure;hold on;
  plot(vessel.center(:,1), vessel.center(:,2), 'r');
  plot(vessel.border(:,1), vessel.border(:,2), 'b');
  plot(vessel.junction.polygon(:,1), vessel.junction.polygon(:,2), 'k')

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

    %{
    figure;hold on;
    scatter(0, 0, 'r');
    scatter(centered(4:end,1), centered(4:end,2), 'b');
    plot([zeros(3,1) centered(1:3,1)].', [zeros(3,1) centered(1:3,2)].', 'g');
    scatter(centered(pos+3,1), centered(pos+3,2), 'k');
    plot([0 centered(away,1)], [0 centered(away,2)], 'k');
    %}

    rads(i) = max(sum(centered(4:end,:).^2, 2));


    %% ERROR here !!!
    vects(:,i) = (-1^(direction(away)==1))*centered(sorting(pos)-1,:) / sqrt(sum(centered(sorting(pos)-1,:).^2));
    junctions(indxs, :) = junctions(polyg, :);
  end

  junct = get_struct('junction');

  junct.polygon = junctions;
  junct.threshold = rads;
  junct.vector = vects;

  return;
end

function centers = polarize_flow(centers)

  [nodes, indxi, indxj] = unique(centers, 'rows');
  goods = ~any(isnan(nodes), 2);

  nodes = nodes(goods, :);
  indxi = indxi(goods, :);

  %figure;hold on
  %plot(centers(:,1), centers(:,2), 'r')
  %scatter(nodes(:,1), nodes(:,2), 'k');

  for i=1:length(indxi)
    connectivity = sum(bsxfun(@eq, indxj(indxi), indxj.'), 2);

    starts = (connectivity == 1);
    nstarts = sum(starts);
    if (nstarts == 0)
      break;
    end

    start_indx = indxi(starts);
    start_indx = start_indx(randi(nstarts, 1));
    nexts = start_indx;

    for j=1:length(indxj)
      %scatter(centers(nexts, 1), centers(nexts, 2), 'g')
      pos = mod(nexts, 3);

      flip = (pos == 2);

      if (any(flip))
        indxs = nexts(flip);

        tmp_pos = centers(indxs, :);
        centers(indxs, :) = centers(indxs-1, :);
        centers(indxs-1, :) = tmp_pos;

        indxi = indxi + ismember(indxi, indxs-1) - ismember(indxi, indxs);

        tmp_indx = indxj(indxs, :);
        indxj(indxs, :) = indxj(indxs-1, :);
        indxj(indxs-1, :) = tmp_indx;

        nexts(flip) = indxs-1;
      end

      tmp_nexts = indxj(nexts+1);
      indxj([nexts nexts+1]) = 0;

      %scatter(nodes(tmp_nexts, 1), nodes(tmp_nexts, 2), 'b')

      nexts = find(ismember(indxj, tmp_nexts));

      if (isempty(nexts))
        break
      end
    end
  end

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
