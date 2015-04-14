function opts = create_vessel(opts)

  b_leachi = get_struct('botrylloides_leachi');
  vessel_props = gmdistribution(b_leachi.vessel_width.mu / opts.pixel_size, b_leachi.vessel_width.sigma / opts.pixel_size, b_leachi.vessel_width.proportions);
  %vessel_widths = random(vessel_props, opts.init_params(1));
  %vessel_centers = randi(prod(opts.image_size), [1 opts.init_params(1)]);

  rim = (opts.image_size - 1) * opts.outside_ridge;
  bounding_box = [1-rim(1) opts.image_size(1)+rim(1) 1-rim(2) opts.image_size(2)+rim(2)];

  max_trials = 200;

  precision = min(opts.image_size)/20;

  vessels = NaN(2,1);
  for i=1:opts.init_params(1)
    %[centery, centerx] = ind2sub(opts.image_size, vessel_centers(i));
    %vessels(i, :) = [rand(1)*(2*pi) vessel_widths(i) centerx centery];

    %get_vessel(bounding_box, vessels(i, 1), vessels(i, 2), vessels(i,3:4), opts.init_params(2));
    for j=1:max_trials
      [new_vessel, new_precision] = get_vessel(bounding_box, vessels, vessel_props, opts.init_params(2));

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

    precision = min(precision, new_precision);

    %p1 = get_bounding_box(vessels(i,:), opts);
    %p2 = get_bounding_box(vessels(i,:)+[pi 0 0 0], opts);
    %[x,y] = polybool('|', p1(1,:), p1(2,:), p2(1,:), p2(2,:));
    %plot(x, y, 'b');
  end

  figure;
  [p,t]=distmesh_poly(vessels.', precision,bounding_box([1 3; 2 4]));

  sort_mesh(p, t);

  %opts.creation_params = vessels(1, 1:2);
  %opts.movement_params = [cos(vessels(1,1)) sin(vessels(1, 1))];

  return;
end

function [vessel, precision] = get_vessel(bounding_box, checks, props, nbranching)

  precision = [];

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
        tmp_len = branches(1, i) * len;
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

      plot(vessel(1,:), vessel(2,:), 'Color', [1 0 0]*(2*i-1)/(2*nbranching))

      vessel = fix_branching(vessel, center, [orig_width curr_width.'], [orig_angle angles]);

      plot(vessel(1,:), vessel(2,:), 'Color', [1 0 0]*(2*i)/(2*nbranching))
      scatter(center(1), center(2), 'r');

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

  precision = min(widths(:)) * 0.75;

  plot(middle(1,:), middle(2,:), 'g')
  plot(vessel(1,:), vessel(2,:), 'b')
  scatter(vessel(1, 1), vessel(2, 1), 'b')
  scatter(vessel(1, 2), vessel(2, 2), 'k')

  %{
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

function vessel = fix_branching(vessel, center, widths, angles)

  ranges = [min(vessel, [], 2) max(vessel, [], 2)];
  vals = [ranges(1, 1):ranges(1, 2)];

  prec = 1024*eps;

  tans = tan(angles);
  lens = NaN(1, 3);
  bkg = NaN(1,3);

  %hfig=figure;hold on;
  %plot(vessel(1,:), vessel(2,:), 'k')
  %scatter(center(1), center(2), 'r');

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

    %plot(vals, vals*tans(i) + bkg(i), 'g');
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

  %scatter(ptsx, ptsy, 'b');
  %scatter(ptsx(bads), ptsy(bads), 'm');
  %scatter(ptsx(in), ptsy(in), 'g');

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

      %if (index == 1)
      %  scatter(new_vessel(1, [1 end-1]), new_vessel(2, [1 end-1]), 'y')
      %else
      %  scatter(new_vessel(1, [index index-1]), new_vessel(2, [index index-1]), 'y')
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
  end

  %plot(new_vessel(1,:), new_vessel(2,:), 'r')

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
