function new_names = register_stack(files, min_frac)

  if nargin<1, files = '*.tif'; min_frac=50; end
  if nargin<2, min_frac=50; end

  new_names = {};
  dir_out = '_registered';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  try
    load_jars();
    turboReg = TurboReg_;
  catch
    error('Registration:register_stack', 'TurboReg is not working properly, please follow the instructions from install_leachi_flow to fix this issue.');
  end

  fprintf(' Registering images pairwise :   .');

  center = ceil(N/2);

  filename = files{center};
  [filepath, fname, fileext] = fileparts(filename);
  new_name = fullfile(out_path, [fname fileext]);
  copyfile(filename, new_name);

  im = imread(filename);
  im = imfillborder(im);

  [hf,wf,c] = size(im);

  is_rgb = (c>1);
  if (is_rgb)
    im = rgb2gray(im);
  end

  max_space = java.lang.Runtime.getRuntime.maxMemory;
  ratio = max_space / (30*2*hf*wf);

  if (ratio < 1)
    im = imresize(im, ratio);
    [h,w,c] = size(im);
  else
    h = hf;
    w = wf;
  end

  xx = repmat(single([1:w]), h, 1);
  yy = repmat(single([1:h]).', 1, w);

  minrad = ceil(max(h,w)/min_frac);
  minsize = (minrad^2);
  hdil = strel('disk', ceil(minrad/20));

  outer = false([h w]);
  outer([1:minrad+1 end-minrad:end],:) = true;
  outer(:, [1:minrad+1 end-minrad:end]) = true;

  [im, pts1] = im2reference(im);

  anchor = fullfile(out_path, 'anchor.tiff');
  target = fullfile(out_path, 'target.tiff');
  source = fullfile(out_path, 'source.tiff');

  imwrite(im, anchor, 'TIFF');
  copyfile(anchor, target);

  command1 =  '-align -file "';
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "'];
  command3 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -rigidBody '];
  commandp = ['%f %f %f %f %f %f %f %f %f %f %f %f -hideOutput'];

  %            num2str(w/2) ' ' num2str(h/2) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
  %            num2str(w/2) ' ' num2str(h/4) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
  %            num2str(w/2) ' ' num2str(3*h/4) ' ' num2str(w/2) ' ' num2str(3*h/4) ...

  new_names = files;
  new_names{center} = new_name;

  for i = center+1:N % loop over images to resize images
    fprintf('\b\b\b%3d', i);

    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    %disp(fname);

    im = imread(filename);
    im = imfillborder(im);
    orig_im = im;

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    [im, pts2] = im2reference(im, pts1);

    imwrite(im, source, 'TIFF');

    command4 = sprintf(commandp, [pts2 pts1].');

    turboReg.run(java.lang.String([command1 source command2 target command3 command4]));
    spts = turboReg.getSourcePoints();
    tpts = turboReg.getTargetPoints();

    H = AffineModel2D(spts(1:3,:), tpts(1:3,:));

    if (ratio < 1)
      H(1:2,3) = H(1:2,3) / ratio;
    end

    im = myimtransform(orig_im, 'affine', H, [wf hf], [0 0]);
    im = imfillborder(im);
    imwrite(im, new_name, 'TIFF');

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    [im, pts1] = im2reference(im);

    imwrite(im, target, 'TIFF');

    new_names{i} = new_name;
  end

  copyfile(anchor, target);

  for i = center-1:-1:1
    fprintf('\b\b\b%3d', i);

    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    %disp(fname);

    im = imread(filename);
    im = imfillborder(im);
    orig_im = im;

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    [im, pts2] = im2reference(im, pts1);

    imwrite(im, source, 'TIFF');

    command4 = sprintf(commandp, [pts2 pts1].');

    turboReg.run(java.lang.String([command1 source command2 target command3 command4]));
    spts = turboReg.getSourcePoints();
    tpts = turboReg.getTargetPoints();

    H = AffineModel2D(spts(1:3,:), tpts(1:3,:));

    if (ratio < 1)
      H(1:2,3) = H(1:2,3) / ratio;
    end

    im = myimtransform(orig_im, 'affine', H, [wf hf], [0 0]);
    im = imfillborder(im);
    imwrite(im, new_name, 'TIFF');

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    [im, pts1] = im2reference(im);

    imwrite(im, target, 'TIFF');

    new_names{i} = new_name;
  end

  delete(target);
  delete(anchor);
  if (N>1)
    delete(source);
  end
  fprintf('\b\b\b\bdone\n');

  return;

  function [img, pts] = im2reference(img, prev_pts)

    img = adjust_tif(img);
    th = graythresh(img(:))*max(img(:));
    mask = (img > th);
    mask = imclose(mask, hdil);
    mask = bwareaopen(mask, minsize);
    mask = imdilate(mask, hdil);
    mask = mask & imfill(~(mask | outer), 'holes');

    if (nargin > 1)
      prev_angle = atan2(diff(prev_pts(1:2,2)), diff(prev_pts(1:2, 1)));
      [pts] = get_landmarks(mask, prev_angle);
    else
      [pts] = get_landmarks(mask);
    end

    if (any(mask(:)))
      img(~mask) = 0;
    end
    img = imnorm(img);

    return;
  end

  function [pts] = get_landmarks(areas, obj_angle)

    x = xx(areas);
    y = yy(areas);

    xbar = mean(x);
    ybar = mean(y);

    x = x - xbar;
    y = y - ybar;

    npts = length(x);

    % Calculate normalized second central moments for the region. 1/12 is
    % the normalized second central moment of a pixel with unit length.
    uxx = sum(x.^2)/npts + 1/12;
    uyy = sum(y.^2)/npts + 1/12;
    uxy = sum(x.*y)/npts;

    % Calculate orientation.
    if (uyy > uxx)
      num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
      den = 2*uxy;
    else
      num = 2*uxy;
      den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
    end
    if (num == 0) && (den == 0)
      angle = 0;
    else
      angle = atan2(num, den);
    end

    if (nargin > 1)
      if (abs(angle - obj_angle) > pi/2)
        angle = angle + pi;
      end
    end

    dist = h/4;
    pts = [cos(angle) sin(angle)]*dist;

    pts = bsxfun(@plus, [-pts; zeros(1,2); pts], [xbar ybar]);

    %figure;imagesc(areas);hold on;scatter(pts(:,1), pts(:,2));

    return;
  end
end

function Ha = AffineModel2D(x1c, x2c)

  angle = (atan2(x1c(3, 2) - x1c(2, 2), x1c(3, 1) - x1c(2, 1)) - ...
          atan2(x2c(3, 2) - x2c(2, 2), x2c(3, 1) - x2c(2, 1)));
  cosa = cos(angle);
  sina = sin(angle);

  Ha = [ cosa sina x2c(1,1) - cosa*x1c(1,1) - sina*x1c(1,2); ...
        -sina cosa x2c(1,2) + sina*x1c(1,1) - cosa*x1c(1,2); ...
          0 0 1];

  return;
end
