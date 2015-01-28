function [all_mvals, all_svals, all_nelems] = mymean(all_vals, dim, indexes)
% MYMEAN computes the mean and standard deviation of the provided data along any
% dimension, ignoring NaNs. It also offers the possibility to provide a grouping
% vector to perform these measurements separately on each subset.
%
%   [MEAN] = MYMEAN(VALS) computes the average of VALS along the first non-singleton
%   dimension.
%
%   [...] = MYMEAN(VALS, DIM) computes the average of VALS along dimension DIM.
%
%   [...] = MYMEAN(VALS, DIM, GROUP) computes the average separately in each GROUP.
%   GROUP should be a vector the same length as the number of elements along DIM,
%   each similar value defining the members of a given GROUP.
%
%   [MEAN, STD] = MYMEAN(...) returns in addition the standard deviation.
%
%   [MEAN, STD, NPTS] = MYMEAN(...) also returns the number of considered points,
%   in each dimension/group depending on the required measurements.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % Initialize some variables
  all_mvals = all_vals;
  all_svals = NaN;
  groups = [];

  % Get the size of the original matrix
  vals_size = size(all_vals);

  % Check if the provided dimension makes sense
  if (nargin == 2)
    if (dim > numel(vals_size))
      return;

    % Put all data in the same group
    else
      indexes = ones(vals_size(dim), 1);
    end

  % Get the default behavior
  elseif (nargin == 1)

    % Find the first non-singleton dimension
    dim = find(vals_size > 1, 1);
    if (isempty(dim))
      dim = 1;
    end

    % And put all data into the same group
    indexes = ones(vals_size(dim), 1);

  % Finally, check again if the dimension makes sense
  elseif (nargin == 3)
    if (dim > numel(vals_size))
      return;
    end
  end

  % In this case, there is no averagin needed, only one data point !
  if (vals_size(dim) == 1)
    all_mvals = all_vals;
    all_svals = zeros(size(all_vals));

    return;

  % Here even worse, no data at all...
  elseif (vals_size(dim) == 0)
    return;
  end

  % The whole trick to be flexible, is that we always work along the first dimension,
  % on a 2D array. We simply permute and reshape our data accordingly !

  % Create the permutation index
  perm_dim = [1:length(vals_size)];
  perm_dim(dim) = 1;
  perm_dim(1) = dim;

  % Permute and reshape in 2D our data such that the measures can be made alond the
  % first dimension !
  all_vals = permute(all_vals, perm_dim);
  all_vals = reshape(all_vals, vals_size(dim), []);

  % Here something went terribly wrong...
  if (numel(indexes) ~= vals_size(dim))
    return;
  end

  % We linearize the indexes as we are now working along the first dimension !
  indexes = indexes(:);

  % Find the individual groups
  groups = unique(indexes).';
  ngroups = length(groups);

  % Now we prepare the output data
  all_mvals = NaN(ngroups, numel(all_vals)/vals_size(dim));
  all_svals = all_mvals;
  all_nelems = all_mvals;

  % Loop over each group
  for i = 1:ngroups

    % The group index
    g = groups(i);

    % Extract to corresponding subset of data
    if (ngroups == 1)
      vals = all_vals;
    else
      vals = all_vals(indexes == g, :);
    end

    % Identify the NaNs
    nans = isnan(vals);

    % Count the number of remaining data
    nelems = sum(~nans, 1);

    % And remove the NaNs
    vals(nans) = 0;

    % Good old-school average !
    mvals = sum(vals, 1) ./ nelems;

    % Replace the NaNs where they belong
    mvals(nelems == 0) = NaN;

    % Store the final output
    all_mvals(i, :) = mvals;

    % Compute the standard deviation only if required
    if (nargout > 1)

      % Same trick here, we do the measurement manually to remove the NaNs.
      svals = bsxfun(@minus, vals, mvals).^2;
      svals(nans) = 0;
      svals = sqrt(sum(svals, 1) ./ (nelems - 1));

      % And replace the NaNs where they belong
      svals(nelems == 0) = NaN;
      svals(nelems == 1) = 0;

      % And store the results
      all_svals(i, :) = svals;

      % All results !
      if (nargout > 2)
        all_nelems(i, :) = nelems;
      end
    end
  end

  % Create the new dimension vector
  vals_size(dim) = ngroups;

  % And reshape plus permute back the data into its original size
  all_mvals = reshape(all_mvals, vals_size(perm_dim));
  all_mvals = ipermute(all_mvals, perm_dim);

  % For all required data !
  if (nargout > 1)
    all_svals = reshape(all_svals, vals_size(perm_dim));
    all_svals = ipermute(all_svals, perm_dim);

    if (nargout > 2)
      all_nelems = reshape(all_nelems, vals_size(perm_dim));
      all_nelems = ipermute(all_nelems, perm_dim);
    end
  end

  return;
end
