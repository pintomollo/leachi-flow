function tmp(num)

  %{
  nmax = 200;
  a = 0.05;

  v = 2 + randn([1, nmax])*3;
  v = [v 20 5 9].';

  v = sort(v);

  counts = [1:length(v)];
  sums = cumsum(v);
  avgs = sums ./ counts.';

  residues = triu(bsxfun(@minus, v, avgs.'));
  S = sqrt(sum(residues.^2, 1) ./ (counts-1));

  prc = tinv(1 - a./(2*counts), counts - 2);
  grubb = sqrt(prc.^2 ./ (counts - 2 + prc.^2));

  thresh = S .* grubb .* (counts-1) ./ sqrt(counts);
  m = max(abs(residues), [], 1);

  m > thresh

  m = mean(v)
  r = v - m;
  s = sqrt(sum(r.^2) / (length(v)-1))
  thresh = s*((length(v)-1)/sqrt(length(v)));
  r > thresh

  return;
  %}

  if (nargin == 0)
    num = 1;
  end

  npts = 600;
  nsamples = 6;

  pos = [1:npts];
  params = abs(bsxfun(@times, (1 + randn(3,num)), [50; 100; 2*pi]));
  params(2,:) = params(2,1);
  %params = abs((1 + randn(3,1)) .* [50; 100; 2*pi]);

  x = repmat(pos, nsamples, 1);

  cosine = zeros(size(x));
  for i=1:num
    cosine = cosine + params(1,i)*cos((i*x/params(2,i))*2*pi + params(3,i));
  end
  %cosine = params(1,1)*cos((x/params(2,1))*2*pi + params(3,1));
  y = cosine + 0.5*min(params(1,:))*randn([nsamples length(x)]);

  %{
  n = estimate_harmonics_number(y);
  figure;plot(y);
  title(n);

  return;

  y = cosine + 2*randn([nsamples length(x)]);
  %[ym, ys] = mymean(y);

  %w = exp(-(y - (ym + ys)).^2 / (ys^2)) + exp(-(y - (ym - ys)).^2 / (ys^2));
  %}

  [period, ampls, phases] = lsqmultiharmonic(x, y);
  nharm = length(ampls);

  cosine = zeros(size(pos));
  for i=1:nharm
    cosine = cosine + ampls(i)*cos((i*pos/period)*2*pi + phases(i));
  end

  figure;scatter(x(:),y(:));
  hold on;plot(x, cosine, 'k');

  [params [ampls period*ones(nharm,1) phases].']

  return;
  keyboard

  p0 = params .* (1 + randn(size(params))*0.2);

  [b,f] = myfit(@err, p0, [0 Inf; 1 Inf; 0 1]);
  keyboard

  return;

  if (nargin==0)
    num=1;
  end

  if (num==1)
    run_orig();
  else
    run_new();
  end

  return;

  function val = err(p, junk)

    val = sum((y - p(1)*sin(((x/p(2)) - p(3))*2*pi) .* w).^2);
    hold off;
    scatter(x,y,'r');
    hold on;
    scatter(x,p(1)*sin(((x/p(2)) - p(3))*2*pi),'b');

    drawnow

    return;
  end
end

function run_orig()

  addpath(absolutepath('./MatPIV 1.7/PTV'));
  addpath(absolutepath('./MatPIV 1.7/masking'));
  addpath(absolutepath('./MatPIV 1.7/src'));
  addpath(absolutepath('./MatPIV 1.7/filters'));
  addpath(absolutepath('./MatPIV 1.7/postprocessing'));

  cd('./MatPIV 1.7/Demo1/')

  keyboard
  [x,y,u,v,snr]=matpiv('im00.bmp','im04.bmp',[64 64;32 32;16 16],...
                       1,0.5,'multinfft',[],'polymask.mat');

  cd('../..');

  rmpath(absolutepath('./MatPIV 1.7/PTV'));
  rmpath(absolutepath('./MatPIV 1.7/masking'));
  rmpath(absolutepath('./MatPIV 1.7/src'));
  rmpath(absolutepath('./MatPIV 1.7/filters'));
  rmpath(absolutepath('./MatPIV 1.7/postprocessing'));

  figure;quiver(x,y,u,v);
  keyboard

  return;
end

function run_new()

  %imgb = imread('./MatPIV 1.7/Demo2/MatPIV000b.bmp');
  %imgc = imread('./MatPIV 1.7/Demo2/MatPIV000c.bmp');

  imgb = imread('./MatPIV 1.7/Demo3/mpim1b.bmp');
  imgc = imread('./MatPIV 1.7/Demo3/mpim1c.bmp');
  maske = load('./MatPIV 1.7/Demo3/polymask.mat');

  maske = ~maske.maske.msk;

  %[x,y,u,v,snr]=matpiv_nfft(imgb,imgc,[64 64;32 32; 16 16],...
  %                     0.5,3,maske,2);
  %[x,y,u,v,snr]=matpiv_nfft(imgb,imgc,[64 64;32 32; 16 16],0.5, 3, maske);
  [x,y,u,v,snr]=matpiv_nfft(imgb,imgc,[128 128; 64 64; 64 64],0.5, 3, maske);

  figure;
  subplot(2,2,1);imagesc(imgb)
  subplot(2,2,2);imagesc(imgc)
  subplot(2,2,3);quiver(x,max(y(:))-y,u,-v);

  keyboard

  return;
end
