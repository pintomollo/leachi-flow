function new_names = register_stack(files, min_frac, reg_type)

  if nargin<1, files = '*.tif'; min_frac=50; reg_type = 'rigidbody';
  elseif nargin<2, min_frac=50; reg_type = 'rigidbody';
  elseif (nargin < 3)
    if (ischar(min_frac))
      reg_type = min_frac;
      min_frac = 50;
    else
      reg_type = 'rigidbody';
    end
  end

  is_affine = (reg_type(1)=='a');
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

  fprintf(' Registering images pairwise :    .');

  center = ceil(N/2);

  filename = files{center};
  [filepath, fname, fileext] = fileparts(filename);
  new_name = fullfile(out_path, [fname fileext]);

  im = imread(filename);
  im = imfillborder(im);

  imwrite(im, new_name, 'TIFF');

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

  [minrad, minsize, hdil] = im2reference([h w], min_frac);
  %minrad = ceil(max(h,w)/min_frac);
  %minsize = (minrad^2);
  %hdil = strel('disk', ceil(minrad/20));

  %outer = false([h w]);
  %outer([1:minrad+1 end-minrad:end],:) = true;
  %outer(:, [1:minrad+1 end-minrad:end]) = true;

  %[im, pts0] = im2reference(im);
  [im, pts0, ptsa] = get_landmarks(im);

  anchor = fullfile(out_path, 'anchor.tiff');
  target = fullfile(out_path, 'target.tiff');
  source = fullfile(out_path, 'source.tiff');

  imwrite(im, anchor, 'TIFF');

  copyfile(anchor, target);
  pts1 = pts0;

  command1 =  '-align -file "';
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "'];
  command3 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1)];

  commandrb = ' -rigidBody ';
  commanda = ' -affine ';

  commandp = ['%f %f %f %f %f %f %f %f %f %f %f %f -hideOutput'];

  commanda = sprintf([commanda commandp], [ptsa ptsa].');

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

    %[im, pts2] = im2reference(im, pts1);
    [im, pts2] = get_landmarks(im, pts1);

    imwrite(im, source, 'TIFF');

    command4 = sprintf(commandp, [pts2 pts1].');

    turboReg.run(java.lang.String([command1 source command2 target command3 commandrb command4]));
    spts = turboReg.getSourcePoints();
    tpts = turboReg.getTargetPoints();

    H = RigidBodyModel2D(spts(1:3,:), tpts(1:3,:));

    if (is_affine)
      im = myimtransform(im, 'affine', H, [w h], [0 0]);
      imwrite(im, source, 'TIFF');

      turboReg.run(java.lang.String([command1 source command2 target command3 commanda]));
      spts = turboReg.getSourcePoints();
      tpts = turboReg.getTargetPoints();

      H = AffineModel2D(spts(1:3,:), tpts(1:3,:)) * H;
    end

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

    %[im, pts1] = im2reference(im);
    [im, pts1] = get_landmarks(im);

    imwrite(im, target, 'TIFF');

    new_names{i} = new_name;
  end

  copyfile(anchor, target);
  pts1 = pts0;

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

    %[im, pts2] = im2reference(im, pts1);
    [im, pts2] = get_landmarks(im, pts1);

    imwrite(im, source, 'TIFF');

    command4 = sprintf(commandp, [pts2 pts1].');

    turboReg.run(java.lang.String([command1 source command2 target command3 commandrb command4]));
    spts = turboReg.getSourcePoints();
    tpts = turboReg.getTargetPoints();

    H = RigidBodyModel2D(spts(1:3,:), tpts(1:3,:));

    if (is_affine)
      im = myimtransform(im, 'affine', H, [w h], [0 0]);
      imwrite(im, source, 'TIFF');

      turboReg.run(java.lang.String([command1 source command2 target command3 commanda]));
      spts = turboReg.getSourcePoints();
      tpts = turboReg.getTargetPoints();

      H = AffineModel2D(spts(1:3,:), tpts(1:3,:)) * H;
    end

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

    %[im, pts1] = im2reference(im);
    [im, pts1] = get_landmarks(im);

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

  function [img, pts, ptsa] = get_landmarks(img, prev_pts)

    [img, areas] = im2reference(img, minsize, hdil);

    if (nargin > 1)
      %if (is_affine)
      %  prev_pts = [mean(prev_pts(1:2,:)); prev_pts(3,:)];
      %end
      prev_angle = atan2(diff(prev_pts(1:2,2)), diff(prev_pts(1:2, 1)));
      [means, angle] = im2moments(areas, xx, yy, prev_angle);
    else
      [means, angle] = im2moments(areas, xx, yy);
    end

    rotmat = [cos(-angle) -sin(-angle); sin(-angle) cos(-angle)];
    %if (is_affine)
    %  pts = [-0.25 -0.25; -0.25 0.25; 0.25 0];
    %else
    %  pts = [-0.25 0; 0 0; 0.25 0];
    %end

    pts = [-0.25 0; 0 0; 0.25 0];
    pts = bsxfun(@plus, (rotmat*bsxfun(@times, pts, [h w]).').', means);

    if (nargout > 2)
      ptsa = [-0.25 -0.25; -0.25 0.25; 0.25 0];
      ptsa = bsxfun(@plus, (rotmat*bsxfun(@times, ptsa, [h w]).').', means);
    end

    %figure;imagesc(areas);hold on;scatter(pts(:,1), pts(:,2));
    %keyboard

    return;
  end
end

function Ha = RigidBodyModel2D(x1c, x2c)

  angle = (atan2(x1c(3, 2) - x1c(2, 2), x1c(3, 1) - x1c(2, 1)) - ...
          atan2(x2c(3, 2) - x2c(2, 2), x2c(3, 1) - x2c(2, 1)));
  cosa = cos(angle);
  sina = sin(angle);

  Ha = [ cosa sina x2c(1,1) - cosa*x1c(1,1) - sina*x1c(1,2); ...
        -sina cosa x2c(1,2) + sina*x1c(1,1) - cosa*x1c(1,2); ...
          0 0 1];

  return;
end

function Ha = AffineModel2D(x1c, x2c)
      
  % A = [x1c(1,1) x1c(2,1) 1 0 0 0
  %     0 0 0 x1c(1,1) x1c(2,1) 1
  %     x1c(1,2) x1c(2,2) 1 0 0 0
  %     0 0 0 x1c(1,2) x1c(2,2) 1
  %     x1c(1,3) x1c(2,3) 1 0 0 0
  %     0 0 0 x1c(1,3) x1c(2,3) 1];
  A = [0 0 0 0 0 0];
  for i=1:3
      A = [A
          x1c(i,:) 1 0 0 0
          0 0 0 x1c(i,:) 1];
  end
  A = A(2:end,:);
      
  % b = [x2c(1,1)
  %     x2c(2,1)
  %     x2c(1,2)
  %     x2c(2,2)
  %     x2c(1,3)
  %     x2c(2,3)];
  b = 0;
  for i=1:3
      b = [b
          x2c(i,:).'];
  end
  b = b(2:end,:);

  %X è il vettore delle incognite 
  %X = [a11
  %     a12
  %     a13
  %     a21
  %     a22
  %     a23];
  X = A\b;

  Ha = [X(1) X(2) X(3)
      X(4) X(5) X(6)
      0 0 1];

  return;
end
