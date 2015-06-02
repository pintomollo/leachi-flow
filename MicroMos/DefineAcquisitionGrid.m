function grid = DefineAcquisitionGrid(parameters)

  nelems = length(parameters.ImageIndexs);
  sizes = factorize_grid(nelems);

  if (nelems == 1)
    trivial
  end

  props = NaN(nelems, 3);

  for i=1:nelems
    strnum = sprintf(parameters.NumberCharactersNumber,parameters.ImageIndexs(i));
    if strcmp(parameters.ImageFormat, '.mat')
        img = load(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum parameters.ImageFormat]));
        img = cell2mat(struct2cell(img));
    else
        img = imread(fullfile(parameters.ImageFolder, [parameters.ImageBaseName strnum parameters.ImageFormat]));
    end
    img = img(1:parameters.ScaleFactor:end,1:parameters.ScaleFactor:end,:);

    if (size(img, 3) > 1)
      img = rgb2gray(img);
    end

    img = double(img);
    edges = imadm(img, false);

    [mval, stdval] = mymean(edges(:));
    props(i, 1) = mval + stdval;
  end

  indexes = [1:nelems];

  best = 1;
  for i=size(sizes, 1):-1:2
    curr_indx = reshape(indexes, sizes(i,[2 1]));

    if (parameters.GridMode==0)
      curr_indx = permute(curr_indx, [2 1]);
    end

    data = props(curr_indx);

    dbl = data(1:end-1,2:end) + data(2:end, 2:end) + data(1:end-1, 1:end-1);
    dbr = data(1:end-1,1:end-1) + data(2:end, 1:end-1) + data(1:end-1, 2:end);
    dul = data(2:end,2:end) + data(1:end-1, 2:end) + data(2:end, 1:end-1);
    dur = data(2:end,1:end-1) + data(1:end-1, 1:end-1) + data(2:end, 2:end);

    keyboard
  end

  return;
end
