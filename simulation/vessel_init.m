function cells = vessel_init(simul, opts)

  density = simul.cell_density * opts.pixel_size^2;

  b_leachi = get_struct('botrylloides_leachi');
  cell_props = gmdistribution(b_leachi.blood_cell.mu, b_leachi.blood_cell.sigma, b_leachi.blood_cell.proportions);

  if (isstruct(simul.creation_params))
    mapping = simul.creation_params.mesh;

    curr_size = mapping.bounding_box;
    curr_n = ceil(density * prod(curr_size));

    pos = bsxfun(@times, rand(curr_n, simul.ndims), curr_size);

    order = mapping.sorted(:,1);
    starts = [0; order(1:end-1,1)];

    props = mapping.proportions;
    regis = mapping.registration;

    bin = bsxfun(@gt, pos(:,1), starts.') & bsxfun(@le, pos(:,1), order(:,1).');

    [junk, indxs] = max(bin, [], 2);
    props = props(indxs, :);

    goods = (pos(:,2) < props(:,1));
    pos = (pos(goods, :));
    indxs = indxs(goods);

    params = mapping.sorted(:,2:end);
    params = params(indxs, :);

    regis = regis(indxs, :);
    props = props(goods, :);

    rel_x = pos(:,1) - starts(indxs);
    rel_y = pos(:,2);

    rel_y1 = params(:,1) - params(:,3) .* (rel_x - params(:,2));
    rel_y2 = params(:,1) - params(:,4) .* (rel_x - params(:,2));

    uppers = (rel_y1 < rel_y);
    lowers = (rel_y2 > rel_y);
    middles = ~(uppers | lowers);

    valids = ~xor(middles, any(regis(:,3:4), 2));

    uppers = uppers & valids;
    lowers = lowers & valids;
    middles = middles & valids;

    middles = middles(:,[1 1]) & regis(:,3:4);

    rel_x(uppers) = rel_x(uppers) + props(uppers,3) - props(uppers,2);
    rel_y(uppers) = props(uppers,1) - rel_y(uppers);

    rel_x(lowers) = rel_x(lowers) + props(lowers,3);

    rel_x(middles(:,1)) = rel_x(middles(:,1)) - props(middles(:,1),4);
    rel_x(middles(:,2)) = props(middles(:,2),2) - rel_x(middles(:,2));

    regis = regis(valids,:);
    origin = mapping.nodes(regis(:,1),:);

    cangle = cos(regis(:,2));
    sangle = sin(regis(:,2));

    rel_x = rel_x(valids);
    rel_y = rel_y(valids);

    x = rel_x .* cangle - rel_y .* sangle;
    y = rel_x .* sangle + rel_y .* cangle;

    pos = [x y] + origin;
    curr_n = size(pos, 1);

    %{
    %%% Display only

    trig = simul.creation_params.mesh.nodes(simul.creation_params.mesh.edges(1,[1:end 1]), :);
    new_trig = bsxfun(@minus, trig, simul.creation_params.mesh.nodes(mapping.registration(1,1),:));

    rot_angle = mapping.registration(1, 2);
    rot = [cos(-rot_angle) -sin(-rot_angle); sin(-rot_angle) cos(-rot_angle)];
    
    new_trig = (rot * new_trig.').';
    ids = (indxs(valids)==1);

    figure;hold on;
    scatter(pos(:,1), pos(:,2));
    scatter(pos(ids,1), pos(ids,2), 'm');
    plot(simul.creation_params.border(:,1), simul.creation_params.border(:,2), 'k');
    plot(simul.creation_params.border(:,1), simul.creation_params.border(:,2), 'k');

    scatter(simul.creation_params.mesh.nodes(:,1), simul.creation_params.mesh.nodes(:,2), 'k');
    plot(trig(:,1), trig(:,2), 'r');
    plot(new_trig(:,1), new_trig(:,2), 'm');
    scatter(rel_x(ids,1), rel_y(ids,1), 'r');

    keyboard
    %}

    %{
    img_size = simul.image_size - 1;
    center = img_size / 2;
    angle = simul.creation_params(1);

    %density = simul.n_cells / (min(img_size)*simul.creation_params(2));
    curr_size = [sqrt(sum(img_size.^2)) * (1 + 2*simul.outside_ridge) simul.creation_params(2)];
    curr_n = ceil(density * prod(curr_size));

    pos = bsxfun(@times, rand(curr_n, simul.ndims)-0.5, curr_size);
    rot = [pos(:,1)*cos(angle) - pos(:,2)*sin(angle), ...
           pos(:,1)*sin(angle) + pos(:,2)*cos(angle)];

    pos = bsxfun(@plus, rot, center);

    cells = [pos random(cell_props, curr_n)];
    cells = cells(~any(cells(:,3:4)<=0, 2), :);
    %}
  else
    curr_size = (simul.image_size-1) * (1 + 2*simul.outside_ridge);
    curr_n = ceil(density * prod(curr_size));

    pos = bsxfun(@minus, bsxfun(@times, rand(curr_n, simul.ndims), curr_size) + 1, ...
                         (simul.image_size-1)*simul.outside_ridge);
  end

  cells = [pos random(cell_props, curr_n)];
  cells = cells(~any(cells(:,3:4)<=0, 2), :);
  cells(:,3) = cells(:,3) / opts.pixel_size;

  %{
  show_vessels(simul)
  hold on;
  scatter(cells(:,1), cells(:,2), 'k')

  keyboard
  %}

  return;
end
