function new_stack = write_stack(files)

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% NEED TO ALSO ALLOW WRITING SPARSE STACKS
  %%% SO NEED TO CHECK IF RGB FILES OR STACK AND
  %%% IF SPARSE
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%


  if nargin<1, files = '*.tif'; end

  new_stack = '';
  dir_out = '_split'; % '../resized';

  [files, out_path] = get_filenames(files, dir_out);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  disp('Reconstructing stack...');

  filename = files{1};
  info = imfinfo(filename);
  im = imread(filename, 'Info', info);

  [h,w,c] = size(im);

  %% http://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-data-as-multi-frame-tiff-format
  type = class(im);
  if (type(1)=='u')
    byte = Tiff.SampleFormat.UInt;
  else
    byte = Tiff.SampleFormat.Int;
  end
  tagstruct = struct('ImageLength', h, ...
                     'ImageWidth', w, ...
                     'Photometric', Tiff.Photometric.MinIsBlack, ...
                     'BitsPerSample', info.BitsPerSample(1), ...
                     'SamplesPerPixel', 1, ...
                     'RowsPerStrip', h, ...
                     'SampleFormat', byte, ...
                     'Compression', Tiff.Compression.None, ...
                     'PlanarConfiguration', Tiff.PlanarConfiguration.(info.PlanarConfiguration));

  new_stack = cell(c, 2);
  for i = 1:c
    new_stack{i,1} = fullfile(out_path, ['stack_' num2str(i) '.tif']);
    tiffobj = Tiff(new_stack{i,1}, 'a');

    tiffobj.setTag(tagstruct);
    tiffobj.write(im(:,:,i));

    new_stack{i,2} = tiffobj;
  end
  fprintf('.');

  try
    for n = 2:N % loop over images to resize images
      filename = files{n};
      im = imread(filename);

      for i = 1:c
        tiffobj = new_stack{i,2};

        tiffobj.writeDirectory();
        tiffobj.setTag(tagstruct);
        tiffobj.write(im(:,:,i));
      end

      fprintf('.');
    end
  catch exception
    for i = 1:c
      new_stack{i, 2}.close();
    end

    throw(exception)
  end

  for i = 1:c
    new_stack{i, 2}.close();
  end

  fprintf(' done !\n');

  new_stack = new_stack(:,1);

  return;
end
