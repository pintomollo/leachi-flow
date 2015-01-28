function [img, noise] = imdenoise(img, rm_bkg, func, varargin)
% IMDENOISE removes noise in a given image using different filtering functions.
%
%   IMG = IMDENOISE(IMG) denoises IMG using a median filter of size 3x3. IMG can be
%   a stack of images, which will then be filtered separately.
%
%   IMG = IMDENOISE(IMG, RM_BKG) if RM_BKG, removes in addition the background signal
%   as estimated (see estimate_noise.m). Default is RM_BKG=false.
%
%   IMG = IMDENOISE(..., FUNC, ARGS) denoises IMG using the function handler FUNC and
%   the corresponding arguments ARGS. Note that no check is performed on ARGS. Note
%   also that for FUNC=@gaussian_mex without ARGS, a default sigma of 0.6 is used [1];
%   for FUNC=@nl_means, an additional arguments corresponding to the estimated level
%   of white noise (i.e. T, see nl_means.m) is passed to FUNC.
%
%   [IMG, NOISE] = IMDENOISE(...) returns in addition the estimated amount of noise
%   present in the image (see estimate_noise.m).
%
% References:
%   [1] Jaensch, S. et al. Automated tracking and analysis of centrosomes in early
%       Caenorhabditis elegans embryos. Bioinformatics. (2010)
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.06.14

  % Default values and input checks
  if (nargin == 0)
    disp('Error: no image provided !')
    img = NaN;
    noise = NaN(0,4);
  elseif (nargin == 1)
    rm_bkg = false;
    func = @median_mex;
  elseif (nargin == 2)
    if islogical(rm_bkg)
      func = @median_mex;
    else
      func = rm_bkg;
      rm_bkg = false;
    end
  end

  % Copy them to alter the list as needed later
  args = varargin;

  % In this particular case, no rm_bkg was provided, which then means that func is
  % part of the arguments
  if (isa(rm_bkg, 'function_handle'))
    args = [{func}, args];
    func = rm_bkg;
    rm_bkg = false;
  end

  % Remove empty cells
  args = args(~cellfun('isempty', args));
  args = args(cellfun(@(x)(isfinite(x) && x>=0), args));

  % Get the image type and convert to double
  img_class = class(img);
  img = double(img);

  % Estimate the noise level in the current image
  noise = estimate_noise(img);

  % Loop over the stack size
  for i = 1:size(img,3)

    % Few educated guesses
    switch func2str(func)

      % Here we follow [1] in case no sigma is provided
      case 'gaussian_mex'
        if (length(args)==0)
          args = {0.6};
        end

      % Here we match the sigma of filtering with the sigma of noise (see nl_means.m)
      case 'nl_means'
        args = [{noise(i,2)}, args];
    end

    % Filter the image
    tmp_img = func(img(:,:,i), args{:});

    % Remove the background if asked
    if (rm_bkg)
      tmp_img = tmp_img - noise(i,1);
      tmp_img(tmp_img < 0) = 0;

      % We thus need to reset the estimated level of background signal
      noise(i,1) = 0;
    end

    % Store the filtered image
    img(:,:,i) = tmp_img;
  end

  % Set the type of the image back to its original one
  img = cast(img, img_class);

  return;
end
