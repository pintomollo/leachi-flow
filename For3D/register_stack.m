function new_names = register_stack(files)

  if nargin<1, files = '*.tif'; end

  new_names = {};

  if (~iscell(files))
    [filepath, filename, fileext] = fileparts(files);
    ls = dir(files);
    ls = clean_dir(ls);

    N = length(ls);

    if (N == 0)
      ls = dir(fullfile(files, '*.tif'));
      ls = clean_dir(ls);
      N = length(ls);
    end

    files = cell([N 1]);

    for i = 1:N
      files{i} = fullfile(filepath, ls(i).name);
    end
  end

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Registering stack...');

  try
    load_jars();
    turboReg = TurboReg_;
  catch
    error('Registration:register_stack', 'TurboReg is not working properly, please follow the instructions from install_leachi_flow to fix this issue.');
  end

  dir_out = '_registered';

  center = ceil(N/2);

  filename = files{center};
  im = imread(filename);
  im = imfillborder(im);

  [filepath, fname, fileext] = fileparts(filename);
  [shorter_path, prev_dir] = fileparts(filepath);

  if (prev_dir(1) == '_')
    out_path = fullfile(shorter_path, dir_out);
  else
    out_path = fullfile(filepath, dir_out);
  end

  if ~isdir(out_path)
    mkdir(out_path);
  end

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

  target = fullfile(out_path, 'target.tiff');
  source = fullfile(out_path, 'source.tiff');

  imwrite(im, target, 'TIFF');

  command1 = '-align -file "';
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "' target ...
              '" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -rigidBody ' ...
              num2str(w/2) ' ' num2str(h/2) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
              num2str(w/2) ' ' num2str(h/4) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
              num2str(w/2) ' ' num2str(3*h/4) ' ' num2str(w/2) ' ' num2str(3*h/4) ...
              ' -hideOutput'];

  new_names = files;

  for i = 1:N % loop over images to resize images
    filename = files{i};
    [filepath, fname, fileext] = fileparts(filename);
    new_name = fullfile(out_path, [fname fileext]);

    disp(fname);

    if (i==center)
      copyfile(filename, new_name);
    else

      im = imread(filename);
      im = imfillborder(im);
      orig_im = im;

      if (is_rgb)
        im = rgb2gray(im);
      end

      if (ratio < 1)
        im = imresize(im, ratio);
      end

      imwrite(im, source, 'TIFF');

      turboReg.run(java.lang.String([command1 source command2]));
      spts = turboReg.getSourcePoints();
      tpts = turboReg.getTargetPoints();

      H = AffineModel2D(spts(1:3,:), tpts(1:3,:));

      if (ratio < 1)
        H(1:2,3) = H(1:2,3) / ratio;
      end

      orig_im = myimtransform(orig_im, 'affine', H, [wf hf], [0 0]);
      imwrite(imfillborder(orig_im), new_name, 'TIFF');
    end
    new_names{i} = new_name;
  end

  delete(target);
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
