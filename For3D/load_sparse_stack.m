function [all_mats, threshs] = load_sparse_stack(files, threshs)

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

  return;
end

function [idx] = mygraythresh(counts)

  % Variables names are chosen to be similar to the formulas in
  % the Otsu paper.
  p = counts / sum(counts);
  omega = cumsum(p);
  mu = cumsum(p .* (1:length(counts))');
  mu_t = mu(end);

  sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

  % Find the location of the maximum value of sigma_b_squared.
  % The maximum may extend over several bins, so average together the
  % locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
  % then return 0.
  maxval = max(sigma_b_squared);
  isfinite_maxval = isfinite(maxval);
  if isfinite_maxval
    idx = mean(find(sigma_b_squared == maxval)) - 1;
  else
    idx = 0;
  end

  return;
end
