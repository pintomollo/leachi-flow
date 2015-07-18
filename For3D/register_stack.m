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

  myMiji(false);
  turboReg = TurboReg_;

  dir_out = 'registered'; % '../resized';

  filename = files{1};
  im = imread(filename);

  [filepath, fname, fileext] = fileparts(filename);
  out_path = fullfile(filepath, dir_out);
  if ~isdir(out_path)
    mkdir(out_path);
  end

  [h,w,c] = size(im);

  if (c > 1)
    im = rgb2gray(im);
  end


  %%%%%%%% JAVA HEAP SPACE MAX SIZE, RESIZE IM ACCORDINGLY

  source = fullfile(out_path, 'source.tiff');
  imwrite(im, source, 'TIFF');

  command1 = ['-align -file "' source '" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -file "'];
  command2 = ['" 0 0 ' num2str(w-1) ' ' num2str(h-1) ' -rigidBody ' ...
              num2str(w/2) ' ' num2str(h/2) ' ' num2str(w/2) ' ' num2str(h/2) ' ' ...
              num2str(w/2) ' ' num2str(h/4) ' ' num2str(w/2) ' ' num2str(h/4) ' ' ...
              num2str(w/2) ' ' num2str(3*h/4) ' ' num2str(w/2) ' ' num2str(3*h/4) ...
              ' -hideOutput'];

  for i = 2:N % loop over images to resize images
    filename = files{i};
    im = imread(filename);

    [filepath, fname, fileext] = fileparts(filename);
    out_path = fullfile(filepath, dir_out);
    if ~isdir(out_path)
      mkdir(out_path);
    end

    [h,w,c] = size(im);

    if (c > 1)
      im = rgb2gray(im);
    end

    target = fullfile(out_path, 'target.tiff');
    imwrite(im, target, 'TIFF');

    keyboard
    turboReg.run(java.lang.String([command1 target command2]));
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
