function new_names = register_stack(files, min_frac)

  if nargin<1, files = '*.tif'; min_frac=100; end
  if nargin<2, min_frac=100; end

  new_names = {};
  dir_out = '_registered';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Registering stack...');

  try
    load_jars();
    turboReg = TurboReg_;
  catch
    error('Registration:register_stack', 'TurboReg is not working properly, please follow the instructions from install_leachi_flow to fix this issue.');
  end

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

  minsize = ceil(max(h,w)/min_frac)^2;
  hdil = strel('disk', 5);

  th = graythresh(im(:))*max(im(:));
  mask = (im > th);
  mask = bwareaopen(mask, minsize);
  mask = imdilate(mask, hdil);

  if (any(mask(:)))
    im(~mask) = 0;
  end

  anchor = fullfile(out_path, 'anchor.tiff');
  target = fullfile(out_path, 'target.tiff');
  source = fullfile(out_path, 'source.tiff');

  imwrite(im, anchor, 'TIFF');
  copyfile(anchor, target);

  command1 =  '-align -file "';
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "'];
  command3 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -rigidBody ' ...
              num2str(w/2) ' ' num2str(h/2) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
              num2str(w/2) ' ' num2str(h/4) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
              num2str(w/2) ' ' num2str(3*h/4) ' ' num2str(w/2) ' ' num2str(3*h/4) ...
              ' -hideOutput'];

  new_names = files;
  new_names{center} = new_name;

  for i = center+1:N % loop over images to resize images
    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    disp(fname);

    im = imread(filename);
    im = imfillborder(im);
    orig_im = im;

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    th = graythresh(im(:))*max(im(:));
    mask = (im > th);
    mask = bwareaopen(mask, minsize);
    mask = imdilate(mask, hdil);

    if (any(mask(:)))
      im(~mask) = 0;
    end

    imwrite(im, source, 'TIFF');

    turboReg.run(java.lang.String([command1 source command2 target command3]));
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

    imwrite(im, target, 'TIFF');

    new_names{i} = new_name;
  end

  copyfile(anchor, target);

  for i = center-1:-1:1
    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    disp(fname);

    im = imread(filename);
    im = imfillborder(im);
    orig_im = im;

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    th = graythresh(im(:))*max(im(:));
    mask = (im > th);
    mask = bwareaopen(mask, minsize);
    mask = imdilate(mask, hdil);

    if (any(mask(:)))
      im(~mask) = 0;
    end

    imwrite(im, source, 'TIFF');

    turboReg.run(java.lang.String([command1 source command2 target command3]));
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

    imwrite(im, target, 'TIFF');

    new_names{i} = new_name;
  end

  delete(target);
  delete(anchor);
  if (N>1)
    delete(source);
  end
  disp('Done !');

  return;
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
