function [nharm, freq] = estimate_harmonics_number(signal, alpha)
% Inspired by Zhang, J. Q. (2003). An eigenvalue residuum-based criterion for detection of the number of sinusoids in white Gaussian noise. Digital Signal Processing: A Review Journal, 13(2), 275–283. doi:10.1016/S1051-2004(02)00029-5
% but rather on the FFT coefficients, and using a Grubb's test instead

  if (nargin < 2)
    alpha = 0.05;
  end

  Fc=fft(signal);
  F=abs(Fc(1:round(end/2)));
  [lambdas, pos] = local_extrema(F(:), 1);

  pos = pos - 1;
  goods = (pos>0);
  lambdas = lambdas(goods);
  pos = pos(goods);

  [lambdas, indxs] = sort(lambdas, 'ascend');

  counts = [1:length(lambdas)];
  sums = cumsum(lambdas);
  avgs = sums ./ counts.';

  residues = triu(bsxfun(@minus, lambdas, avgs.'));
  S = sqrt(sum(residues.^2, 1) ./ (counts-1));

  prc = tinv(1 - alpha./(2*counts), counts - 2);
  grubb = sqrt(prc.^2 ./ (counts - 2 + prc.^2));

  thresh = S .* grubb .* (counts-1) ./ sqrt(counts);

  goods = (max(residues, [], 1) > thresh);

  indx = find(~goods, 1, 'last');
  nharm = length(goods) - indx;

  if (nargout > 1)
    freq = pos(indxs(indx:end))/length(signal);
  end

  return;
end
