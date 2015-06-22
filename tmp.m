function tmp(num)

  x = [1:600];

  params = rand(3,1) .* [50; 30; 1]
  y = params(1)*sin(((x/params(2)) - params(3))*2*pi) + randn([1 length(x)]);
  %[ym, ys] = mymean(y);

  %w = exp(-(y - (ym + ys)).^2 / (ys^2)) + exp(-(y - (ym - ys)).^2 / (ys^2));

  figure;scatter(x,y);

  lsqmultiharmonic(x, y, 2);

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
