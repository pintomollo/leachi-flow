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
