function register_stack(files)

  if nargin<1, files = '*.tif'; end

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

  try
    myMiji(false);
    turboReg = TurboReg_;
  catch
    %%% DISPLAY ERROR
  end

  dir_out = 'registered'; % '../resized';

  filename = files{1};
  im = imread(filename);
  im = imfillborder(im);

  [filepath, fname, fileext] = fileparts(filename);
  out_path = fullfile(filepath, dir_out);
  if ~isdir(out_path)
    mkdir(out_path);
  end

  [h,w,c] = size(im);

  is_rgb = (c>1);
  if (is_rgb)
    im = rgb2gray(im);
  end


  %%%%%%%% JAVA HEAP SPACE MAX SIZE, RESIZE IM ACCORDINGLY
  max_space = java.lang.Runtime.getRuntime.maxMemory;
  ratio = max_space / (30*2*h*w);

  if (ratio < 1)
    im = imresize(im, ratio);
    [h,w,c] = size(im);
  end

  source = fullfile(out_path, 'source.tiff');
  imwrite(im, source, 'TIFF');

  command1 = ['-align -file "' source '" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "'];
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -rigidBody ' ...
              num2str(w/2) ' ' num2str(h/2) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
              num2str(w/2) ' ' num2str(h/4) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
              num2str(w/2) ' ' num2str(3*h/4) ' ' num2str(w/2) ' ' num2str(3*h/4) ...
              ' -hideOutput'];

  command3 = ['-transform -file "' source '" ' num2str(w) ' ' num2str(h) ' -rigidBody '];

  for i = 2:N % loop over images to resize images
    filename = files{i};
    im = imread(filename);
    im = imfillborder(im);

    [filepath, fname, fileext] = fileparts(filename);
    out_path = fullfile(filepath, dir_out);
    if ~isdir(out_path)
      mkdir(out_path);
    end

    if (is_rgb)
      im = rgb2gray(im);
    end

    if (ratio < 1)
      im = imresize(im, ratio);
    end

    target = fullfile(out_path, 'target.tiff');
    imwrite(im, target, 'TIFF');

    turboReg.run(java.lang.String([command1 target command2]));
    spts = turboReg.getSourcePoints();
    tpts = turboReg.getTargetPoints();

    angle = -(atan2(tpts(3,1) - tpts(2,1), tpts(3, 2) - tpts(2, 2)) - ...
            atan2(spts(3, 1) - spts(2, 1), spts(3, 2) - spts(2, 2)));
    cosa = cos(angle);
    sina = sin(angle);

    trans_mat = [ cosa sina spts(1,1) - cosa*tpts(1,1) + sina*tpts(1,2); ...
                 -sina cosa spts(1,2) - sina*tpts(1,1) - cosa*tpts(1,2); ...
                 0 0 1];

    im2 = imtransform(im, maketform('affine', trans_mat.'), 'Size', [h w]);

    modified = fullfile(out_path, 'modified.tiff');
    imwrite(im2, modified, 'TIFF');

    command4 = [num2str(spts(1,1)) ' ' num2str(spts(1,2)) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
                num2str(spts(2,1)) ' ' num2str(spts(2,2)) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
                num2str(spts(3,1)) ' ' num2str(spts(3,2)) ' ' num2str(w/2) ' ' num2str(3*h/4) ...
                ' -hideOutput'];
    turboReg.run(java.lang.String([command3 command4]));
    im3 = turboReg.getTransformedImage();
    im3 = reshape(im3.getProcessor().getPixels(), [w h]).';

    modified = fullfile(out_path, 'turbo_wrap.tiff');
    imwrite(all2uint16(im3), modified, 'TIFF');

    keyboard
  end

  % Not fitting in memory, need to look in StackReg.java, see how it's done and call image by image to do it ourselves using directly pairs of images and TurboReg.

  % run("Image Sequence...", "open=/Users/blanchou/Documents/MATLAB/Movies/Sectioning/Bleachi_Colony/resized");

  

  return;
end

function list = clean_dir(list)

  list = list(~[list.isdir]);
  for i=length(list):-1:1
    if (list(i).name(1)=='.')
      list(i) = [];
    end
  end

  return;
end
