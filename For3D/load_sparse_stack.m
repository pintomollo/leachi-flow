function [all_mats, threshs, type] = load_sparse_stack(files, threshs)

  if nargin<1
    files = '*.tif';
    threshs = [];
  elseif nargin < 2
    threshs = [];
  end

  [files] = get_filenames(files);

  N = length(files);
  if (N == 0), disp('nada??'), return, end

  if (isempty(threshs))
    threshs = NaN([1, N]);
  end

  nbins = 256;
  all_mats = cell([N, 1]);

  for i=1:N
    [nplanes, imgsize] = size_data(files{i});

    smat = sparse([],[],[],prod(imgsize),nplanes,0);

    if (~isfinite(threshs(i)))
      counts = zeros(nbins, 1);
      tiffobj = Tiff(files{i}, 'r');
      for n=1:nplanes-1
        img = tiffobj.read();
        tiffobj.nextDirectory();
        h = imhist(img(:), nbins);
        counts = counts + h;
      end
      img = tiffobj.read();
      tiffobj.close();

      h = imhist(img(:), nbins);
      counts = counts + h;

      threshs(i) = mygraythresh(counts);
    end

    tiffobj = Tiff(files{i}, 'r');
    for n=1:nplanes-1
      img = tiffobj.read();
      tiffobj.nextDirectory();

      img(img<=threshs(i)) = 0;
      smat(:,n) = sparse(double(img(:)));
    end
    img = tiffobj.read();
    tiffobj.close();

    img(img<=threshs(i)) = 0;
    smat(:,nplanes) = sparse(double(img(:)));

    smat = ndSparse(smat, [imgsize nplanes]);
    all_mats{i} = smat;
  end

  if (N==1)
    all_mats = all_mats{1};
  end

  if (nargout>2)
    type = class(img);
  end

  return;
end
