function [gammas] = estimate_noise(img, filter_type, block_size, mfilter, precision)
% ESTIMATE_NOISE returns an estimation of the noise present in the image as:
% [BKG STD LAMBDA MU] where BKG and STD are the mean and standard deviation of
% the uniform white noise (i.e. the thermal noise from the camera), LAMBDA is the
% variance of the linear noise (with respect to intensity, i.e. the Poisson shot noise
% from the emitted photons) and MU is the variance of the multiplicative noise
% (i.e. the quadratic noise from an EMCCD amplification). The algorithm is implemented
% based on [1], which itself is based on [2]. In addition, we have improved [2] by
% replacing their edge detector with IMADM (help imadm).
%
%   NOISES = ESTIMATE_NOISE(IMG) estimates the four sources of NOISES present in IMG.
%
%   NOISES = ESTIMATE_NOISE(STACK) estimates the noises in each plane separately,
%   returning a Nx4 matrix of NOISES.
%
%   NOISES = ESTIMATE_NOISE(..., FILTER) allows to specify the edge detector utilized.
%   Implemented options are: 'adm', 'published' or a custom filtering matrix. Default
%   is 'adm'.
%
%   NOISES = ESTIMATE_NOISE(..., FILTER, BLOCK_SIZE, RADIUS, PRECISION) allows to
%   specify the parameters of the algorithm: BLOCK_SIZE is the size of the tiles the
%   image is partitioned into (see [1]), RADIUS the size of the kernel used by the median
%   filter (as a square) to obtain the "noiseless" image and PRECISION a vector [CI, ER]
%   defining the aimed confidence and percentage of error in the estimation of the
%   standard deviation as defined in [3]. Default values are 21 for BLOCK_SIZE, 3 for
%   NITER and [0.8 0.1] for PRECISION. Provide an empty value to use the default value
%   for any of the parameters.
%
% References:
%   [1] Paul P, Duessmann H, Bernas T, Huber H, Kalamatianos D, Automatic noise
%       quantification for confocal fluorescence microscopy images.
%       Comput Med Imaging Graph 34 (2010)
%   [2] Amer A, Dubois E, Fast and reliable structure-oriented video noise estimation.
%       IEEE Trans Circuits Syst Video Technol 15 (2005)
%   [3] Greenwood J, Sandomire M, Sample size required for estimating the standard
%       deviation as a per cent of its true value. J Am Stat 45 (1950)
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % Input checking
  if(nargin < 2)
    filter_type = '';
  end
  if (nargin < 3)
    block_size = [];
  end
  if (nargin < 4)
    mfilter = [];
  end
  if (nargin < 5)
    precision = [];
  end

  % There is nothing to parse
  if (nargin  == 0 || isempty(img))
    gammas = NaN(1,4);
    return;
  end

  % Default values
  if (isempty(filter_type))
    filter_type = 'adm';
  end
  if (isempty(block_size))
    block_size = 21;
  end
  if (isempty(mfilter))
    mfilter = 3;
  end
  if (numel(precision) < 2)
    precision = [0.8 0.1];
  end

  % In case we have a stack, parse each plane separately
  nplanes = size(img, 3);
  if (nplanes > 1)
    gammas = NaN(nplanes, 4);
    for i=1:nplanes
      gammas(i,:) = estimate_noise(img(:,:,i), filter_type, block_size, mfilter, precision);
    end

    return;
  end

  % Checking if we have a custom filter or not
  if (~ischar(filter_type))
    filter = filter_type;
    filter_type = 'custom';
  end

  % Infer the minimal number of samples required to meet the confidence coefficient
  % when infering the standard deviation
  minsize = estimate_sample_size(precision(1), precision(2));

  % Deduce the number of bins we can use
  nbins = numel(img) / (50*minsize);

  % The number of elements to avoid incomplete blocks on the border of the image
  nelems = block_size^2;

  % Switch filter types
  switch filter_type

    % The default value, fitler using imadm and runs the block processing
    case 'adm'
      nelems = 2*nelems;
      edge_map = imadm(img, 0, false);
      blocks = blockproc(cat(3,img,edge_map), [block_size block_size], @analyze_filt);

    % Similar but using a custom filter
    case 'custom'
      nelems = 2*nelems;
      % We take the absolute value, just in case
      edge_map = abs(imfilter(img, filter, 'symmetric'));
      blocks = blockproc(cat(3,img,edge_map), [block_size block_size], @analyze_filt);

    % The published alternative, it creates a large filter to apply on each
    % block to get a faster estimate of the edge, and thus of the uniformity
    otherwise

      % The index of the middle pixel
      middle = ((block_size-1)/2)+1;

      % Now let's build these filters, see [2], but basically, a vertical, an horizontal
      % two diagonaly and four corners with each their proper weights
      row = -ones(1, block_size);
      hfilter = zeros(block_size);
      hfilter(middle, :) = row;
      vfilter = hfilter.';
      pdfilter = diag(row);
      ndfilter = pdfilter.';
      ldfilter = hfilter;
      ldfilter(middle,middle:end) = 0;
      ldfilter(middle:end, middle) = -1;
      lufilter = hfilter;
      lufilter(middle,middle:end) = 0;
      lufilter(1:middle-1, middle) = -1;
      rdfilter = hfilter;
      rdfilter(middle,1:middle) = 0;
      rdfilter(middle+1:end, middle) = -1;
      rufilter = hfilter;
      rufilter(middle,1:middle) = 0;
      rufilter(1:middle-1, middle) = -1;

      % Finally, concatenate all filters
      filter = [hfilter(:) vfilter(:) pdfilter(:) ndfilter(:) ldfilter(:) lufilter(:) rdfilter(:) rufilter(:)];
      % And adjust their weights
      filter((middle-1)*block_size + middle, :) = (block_size - 1);

      % The runs the block processing, with a slightly different function for speedup
      blocks = blockproc(img, [block_size block_size], @analyze_block);
  end

  % We use the median filter to infer the noiseless image
  noisefree = median_mex(img, mfilter);

  % Extract the noise only
  noisy = img - noisefree;

  % Reshape the output for simplicity
  blocks = reshape(blocks, [], 3);

  % We sort the resulty to find the most homogeneous blocks, with the lowest signal
  blocks = sortrows(blocks);

  % And we define the target standard deviation on the most uniform blocks
  target_var = mean(blocks(1:3, 3));

  % As defined in [1], we keep only the low variance blocks
  goods = (blocks(:,3) <= 3*target_var);
  blocks = blocks(goods, :);

  % And we estimate the background and the standard deviation on the remaining blocks,
  % with a limit on their number to avoid too non-uniform blocks if we can afford it
  nblocks = min(size(blocks, 1), 20);
  gauss_noise = mean(blocks(1:nblocks, 2:3));

  % Now basically, we want to perform a linear regression on the amount of noise
  % in each pixel, based on its noiseless intensity (see [1])
  noisefree = noisefree(:);
  noisy = noisy(:);

  % Build the histogram edges
  edges = ([0:nbins].')*range(noisefree)/nbins + min(noisefree);
  edges(1) = edges(1) - 1e-6;
  edges(end) = edges(end) + 1e-6;

  % Now classify the pixel noiseless intensities, keep the map for identifying the
  % corresponding noises
  [counts, map] = histc(noisefree, edges);

  % Initialize the variables for the variance estimation
  nbins = length(counts);
  vars = NaN(nbins, 1);

  % We now loop over the histogram, keeping only the bins that contain enough data
  % as defined by [3]. For those, we estimate the standard deviation between all
  % pixels which have the similar intensity
  for i=1:nbins
    if (counts(i) > minsize)
      vars(i) = std(noisy(map==i))^2;
    end
  end

  % We keep only the variance we actually estimates, and we get rid of very high
  % intensities as they are usually biased
  goods = isfinite(vars);
  goods(ceil(end/2):end) = false;
  vars = vars(goods);

  % Center the pixel intensity for the upcoming linear regression
  edges = (edges(1:end-1) + edges(2:end)) / 2;
  edges = edges(goods);

  % Correct for the white noise
  edges = edges - gauss_noise(1);
  vars = vars - gauss_noise(2)^2;

  % Perform the actual regression
  X = [edges edges.^2];
  coefs = X \ vars;

  % Finally, we need to check whether some of the coefficients do not make sense
  switch (sum(coefs <= 0))

    % Everything is fine, concatenate and output
    case 0
      gammas = [gauss_noise coefs.'];

    % There was some overfitting, usually because there is no multiplicative noise,
    % we simply refit the correct data alone and set the other one to 0
    case 1
      if (coefs(1) < 0)
        coefs = [edges.^2] \ vars;
        gammas = [gauss_noise 0 coefs];
      else
        coefs = [edges] \ vars;
        gammas = [gauss_noise coefs 0];
      end

      % Just in case something went wrong twice
      gammas(gammas < 0) = 0;

    % No additional data...
    case 2
      gammas = [gauss_noise 0 0];
  end

  return;

  function res = analyze_filt(block_struct)
  % The block-processing function which estimates the homogeneity, the average and
  % the standard deviation in each block.

    % Get the data
    data = block_struct.data(:);

    % Run only if we have a whole block
    if (numel(data) == nelems)

      % We compute the homogeneity based on the edge intensity
      homog = sum(data((end/2)+1:end));
      [means, stds] = mymean(data(1:end/2));

      % And format the output
      res = cat(3, homog, means, stds);

    % Otherwise, we use the highest value as later we'll be looking for low-value
    % blocks !
    else
      res = Inf([1 1 3]);
    end

    return;
  end

  function res = analyze_block(block_struct)
  % Almost the same as analyze_filt, except that here we need to apply the custom
  % filter developed in [2] to the data first.

    % Get the data
    data = block_struct.data(:);

    % Check if it's a full block
    if (numel(data) == nelems)

      % Filter the data and get the other measurements
      homog = sum(abs(sum(bsxfun(@times, data, filter), 1)), 2);
      [means, stds] = mymean(data);

      % The output
      res = cat(3, homog, means, stds);

    % A problem
    else
      res = Inf([1 1 3]);
    end

    return;
  end
end

function sample_size = estimate_sample_size(ci, precision)
% This function estimates the number of samples required to estimate a standard
% deviation which meets a given confidence coefficient and reaches a certain level
% of precision (as a percentage of its value). See [3] for more information.

  % An non-exhaustive list of sample number
  n = 10.^[0:4];

  % The threshold precision we want to obtain
  p_plus = (1+precision)^2;
  p_minus = (1-precision)^2;

  % The precision at the initial sizes
  ps = curr_prob(n);

  % Check if we can meet the criterion (i.e. is there a 0 somewhere ?)
  cross_pos = (abs(diff(sign(ps))) > 0); 

  % If yes, we use fzero to find the exact number
  if (length(cross_pos) > 0 & any(cross_pos))
    indx = find(cross_pos, 1);
    ns = n([indx indx+1]);

    % The sample size !
    N = fzero(@curr_prob, ns);
    sample_size = ceil(N) + 1;

  % Otherwise, we need to give up
  else
    if (all(ps < 0 | ~isfinite(ps)))
      sample_size = Inf;
    else
      sample_size = 1;
    end
  end

  return;

  function prob_precision = curr_prob(dof)
  % This function derives directly from [3], and allows one to directly estimate
  % the precision given a certain number of samples.

    chi_high = dof*p_plus;
    chi_low = dof*p_minus;

    p1 = chi2cdf(chi_high, dof);
    p2 = chi2cdf(chi_low, dof);

    prob_precision = ci - (p1 - p2);

    return;
  end
end
