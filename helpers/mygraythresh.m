function [idx] = mygraythresh(counts, type)

  % In case we are provided directly with a sparse matrix
  if (issparse(counts))

    % We need to know the internal data type
    if (nargin < 2)
      type = 'uint8';
    end

    % Create the bins for the histogram
    maxval = double(intmax(type));
    nbins = 256;
    edges = [0:nbins] * (maxval / (nbins-1));

    % The number of empty elements
    nzeros = numel(counts) - nnz(counts);

    % Counts the elements
    counts = histc(nonzeros(counts(:)), edges);
    counts = counts(1:end-1);
    counts(1) = counts(1) + nzeros;
  end

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
