function [mapped] = LKTracker( img1, img2, pts, shift)
% [ X2 Y2 ] = LKTrackPyr( img1, img2, X1, Y1 )
%	Lucas-Kanade Tracker with pyramid and iteration.
%	img1 and img2 are 2 images of a same scene with little time-lag,
%	X's and Y's are coordinates to be tracked.

% Adapted from http://www.mathworks.com/matlabcentral/fileexchange/30822

  winR = 5;
  th = .5;
  %th = .01;
  maxIter = 20;
  minImgSz = 64; % if pyramid level is too high, corners may be blurred
  maxPyrLevel = floor(log2(min(size(img1))/minImgSz));
  scaling = 1 / 2^(maxPyrLevel);

  [winX,winY] = meshgrid(-winR:winR,-winR:winR);

  img1pyrs = genPyr(img1,maxPyrLevel);
  img2pyrs = genPyr(img2,maxPyrLevel);
  h = fspecial('sobel');

  X1 = pts(:,1);
  Y1 = pts(:,2);

  X2 = X1 - shift(1);
  Y2 = Y1 - shift(2);

  ptNum = size(X1,1);

  X1 = X1 * scaling;
  Y1 = Y1 * scaling;
  X2 = X2 * scaling;
  Y2 = Y2 * scaling;

  updated = false(size(X1));

  %X2 = X1/2^maxPyrLevel;
  %Y2 = Y1/2^maxPyrLevel;

  for level = maxPyrLevel:-1:1
    img1 = img1pyrs{level};
    img2 = img2pyrs{level};
    [M N] = size(img1);

    img1x = imfilter(img1,h','replicate');
    img1y = imfilter(img1,h,'replicate');

    X1 = X1*2;
    Y1 = Y1*2;
    X2 = X2*2;
    Y2 = Y2*2;

    for p = 1:ptNum
      %xt = X1(p);
      %yt = Y1(p);

      %[iX iY oX oY isOut] = genMesh(xt,yt,winR,M,N);
      [oX oY isOut] = tmp_genMesh(X1(p),Y1(p));

      if isOut
        continue;
      end % if the window of a point is out of image bound

      Ix = bilinear_mex(img1x,oX,oY); % notice the order of X and Y
      Iy = bilinear_mex(img1y,oX,oY);
      I1 = bilinear_mex(img1,oX,oY);

      %Ix = interp2(iX,iY,img1x(iY,iX),oX,oY); % notice the order of X and Y
      %Iy = interp2(iX,iY,img1y(iY,iX),oX,oY);
      %I1 = interp2(iX,iY,img1(iY,iX),oX,oY);

      mat = [Ix(:) Iy(:)];
      if (any(isnan(mat(:))) || rank(mat) < 2)
        continue;
      end

      status = false;
      xt = X2(p);
      yt = Y2(p);

      for q = 1:maxIter
        [oX oY isOut] = tmp_genMesh(xt,yt);
        %[iX iY oX oY isOut] = genMesh(xt,yt,winR,M,N);

        if isOut
          break;
        end

        It = bilinear_mex(img2,oX,oY) - I1;

        vel = mat \ It(:);
        xt = xt+vel(1);
        yt = yt+vel(2);

        if max(abs(vel))<th || isnan(vel(1))
          status = true;
          break;
        end
      end

      if (~isnan(xt))
        X2(p) = xt;
        Y2(p) = yt;

        updated(p) = status;
      end
    end
  end

  X2(~updated) = M+1;
  Y2(~updated) = N+1;

  mapped = [X2 Y2];

  return;

  function [coX, coY, cOut] = tmp_genMesh(xt,yt)

    coX = winX + xt;
    coY = winY + yt;

    cOut = (coX(1)<1 || coY(1)<1 || coX(end)>N || coY(end)>M);

    return;
  end
end

function [iX iY oX oY isOut] = genMesh(xt,yt,winR,M,N)

  l = xt-winR;
  t = yt-winR;
  r = xt+winR;
  b = yt+winR;

  [oX,oY] = meshgrid(l:r,t:b);
  fl = floor(l);
  ft = floor(t);
  cr = ceil(r);
  cb = ceil(b);
  iX = fl:cr;
  iY = ft:cb;

  isOut = (fl<1 || ft<1 || cr>N || cb>M);
  % error('out of image bound.');

  return;
end

%function [ pyr ] = genPyr( img, type, level )
function [ pyr ] = genPyr( img, level )
%[ pyr ] = genPyr( img, type, level )
%	generate Gaussian or Laplacian pyramid
%   PYR = GENPYR(A,TYPE,LEVEL) A is the input image, 
%	can be gray or rgb, will be forced to double. 
%	TYPE can be 'gauss' or 'laplace'.
%	PYR is a 1*LEVEL cell array.

  %{
  pyr = cell(1,level);
  pyr{1} = im2double(img);

  for p = 2:level
          pyr{p} = pyr_reduce(pyr{p-1});
  end

  if strcmp(type,'gauss'), return; end

  for p = level-1:-1:1 % adjust the image size
          osz = size(pyr{p+1})*2-1;
          pyr{p} = pyr{p}(1:osz(1),1:osz(2),:);
  end

  for p = 1:level-1
          pyr{p} = pyr{p}-pyr_expand(pyr{p+1});
  end
  %}

  pyr = cell(1,level);
  img = im2double(img);

  pyr{1} = img;

  kernelWidth = 5; % default
  cw = .375; % kernel centre weight, same as MATLAB func impyramid. 0.6 in the Paper
  ker1d = [.25-cw/2 .25 cw .25 .25-cw/2];
  kernel = kron(ker1d,ker1d');

  for p = 2:level
    imgFiltered = imfilter(img,kernel,'replicate','same');
    img = imgFiltered(1:2:end,1:2:end,:);

    pyr{p} = img;
  end

  return;
end

%{
function [ imgout ] = pyr_reduce( img )
  %PYR_REDUCE  Image pyramid reduction
  %   B = PYR_REDUCE( A )  If A is M-by-N, then the size of B 
  %	is ceil(M/2)-by-ceil(N/2). Support gray or rgb image.
  %	B will be transformed to double class.
  %	Results the same w/ MATLAB func impyramid.

  kernelWidth = 5; % default
  cw = .375; % kernel centre weight, same as MATLAB func impyramid. 0.6 in the Paper
  ker1d = [.25-cw/2 .25 cw .25 .25-cw/2];
  kernel = kron(ker1d,ker1d');

  img = im2double(img);
  sz = size(img);
  imgout = [];

  for p = 1:size(img,3)
          img1 = img(:,:,p);
          imgFiltered = imfilter(img1,kernel,'replicate','same');
          imgout(:,:,p) = imgFiltered(1:2:sz(1),1:2:sz(2));
  end

end
%}
