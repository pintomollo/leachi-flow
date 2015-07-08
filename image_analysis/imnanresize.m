function out = imnanresize(img, output_size)

  kernel = @cubic;
  kernel_size = 4;
  [m,n,o] = size(img);
  img_size = [m n o];
  scale = output_size ./ img_size(1:2);
  antialiasing = true;

  % Calculate interpolation weights and indices for each dimension.
  weights = cell(1,2);
  indices = cell(1,2);

  out = img;
  new_size = img_size;
  for k = 1:2
    [weights{k}, indices{k}] = contributions(img_size(k), ...
        output_size(k), scale(k), kernel, ...
        kernel_size, antialiasing);

    new_size(k) = output_size(k);

    in = out;
    out = NaN(new_size);

    for p=1:new_size(3)
      for i=1:new_size(1)
        for j=1:new_size(2)
          if (k==1)
            tmp_vals = in(indices{k}(i,:),j,p);
            tmp_weig = weights{k}(i,:).';

            %out(i,j) = sum(in(indices{k}(i,:),j).*weights{k}(i,:).', k);
          else
            tmp_vals = in(i,indices{k}(j,:),p);
            tmp_weig = weights{k}(j,:);

            %out(i,j) = sum(in(i,indices{k}(j,:)).*weights{k}(j,:), k);
          end
          bads = isnan(tmp_vals);

          if (all(bads))
            out(i,j,p) = NaN;
          else
            tmp_vals(bads) = 0;
            tmp_weig(bads) = 0;
            tmp_weig = tmp_weig / sum(tmp_weig);

            out(i,j,p) = sum(tmp_vals.*tmp_weig, k);
          end
        end
      end
    end
  end

  return;
end

function [weights, indices] = contributions(in_length, out_length, ...
                                            scale, kernel, ...
                                            kernel_width, antialiasing)


  if (scale < 1) && (antialiasing)
      % Use a modified kernel to simultaneously interpolate and
      % antialias.
      h = @(x) scale * kernel(scale * x);
      kernel_width = kernel_width / scale;
  else
      % No antialiasing; use unmodified kernel.
      h = kernel;
  end

  % Output-space coordinates.
  x = (1:out_length)';

  % Input-space coordinates. Calculate the inverse mapping such that 0.5
  % in output space maps to 0.5 in input space, and 0.5+scale in output
  % space maps to 1.5 in input space.
  u = x/scale + 0.5 * (1 - 1/scale);

  % What is the left-most pixel that can be involved in the computation?
  left = floor(u - kernel_width/2);

  % What is the maximum number of pixels that can be involved in the
  % computation?  Note: it's OK to use an extra pixel here; if the
  % corresponding weights are all zero, it will be eliminated at the end
  % of this function.
  P = ceil(kernel_width) + 2;

  % The indices of the input pixels involved in computing the k-th output
  % pixel are in row k of the indices matrix.
  indices = bsxfun(@plus, left, 0:P-1);

  % The weights used to compute the k-th output pixel are in row k of the
  % weights matrix.
  weights = h(bsxfun(@minus, u, indices));

  % Normalize the weights matrix so that each row sums to 1.
  weights = bsxfun(@rdivide, weights, sum(weights, 2));

  % Clamp out-of-range indices; has the effect of replicating end-points.
  indices = min(max(1, indices), in_length);

  % If a column in weights is all zero, get rid of it.
  kill = find(~any(weights, 1));
  if ~isempty(kill)
      weights(:,kill) = [];
      indices(:,kill) = [];
  end
end

function f = cubic(x)
% See Keys, "Cubic Convolution Interpolation for Digital Image
% Processing," IEEE Transactions on Acoustics, Speech, and Signal
% Processing, Vol. ASSP-29, No. 6, December 1981, p. 1155.

  absx = abs(x);
  absx2 = absx.^2;
  absx3 = absx.^3;

  f = (1.5*absx3 - 2.5*absx2 + 1) .* (absx <= 1) + ...
                  (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) .* ...
                  ((1 < absx) & (absx <= 2));

  return;
end
