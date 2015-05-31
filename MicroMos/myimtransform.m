function [img] = myimtransform(img, modello, GLOBAL, XData, YData)

  is_rgb = (size(img, 3) > 1);

  if (nargin < 5)
    XData = modello(1):modello(2);
    YData = GLOBAL(1):GLOBAL(2);

    goodx = (XData >= 0 & XData < size(img,2));
    goody = (YData >= 0 & YData < size(img,1));

    if (is_rgb)
      tmp_img = NaN(length(YData), length(XData), 3);
      tmp_img(goody, goodx, :) = img;
    else
      tmp_img = NaN(length(YData), length(XData));
      tmp_img(goody, goodx) = img;
    end
    img = tmp_img;
  else

    [X,Y] = meshgrid(XData(1):XData(2),YData(1):YData(2));
    indxs = [X(:) Y(:) zeros(numel(X), 1)] + 1;
    %Tinv = inv(GLOBAL');

    if (modello(1)=='a')
      %Tinv(1:end-1,end) = 0;
      %Tinv(end,end) = 1;

      U1 = indxs / GLOBAL';                  % Transform in homogeneous coordinates
      inv_indxs  = U1(:,1:end-1);           % Convert homogeneous coordinates to U

      %T = maketform(modello,GLOBAL');
      %[P2] = T.inverse_fcn([X(:) Y(:)]+1, T);
    else
      U1 = indxs / GLOBAL';                  % Transform in homogeneous coordinates
      inv_indxs  = bsxfun(@rdivide, U1(:,1:end-1), U1(:,end));     % Convert homogeneous coordinates to U
    end

    if (is_rgb)
      img_r = bilinear_mex(double(img(:,:,1)), inv_indxs);
      img_r = reshape(img_r, size(X));
      img_g = bilinear_mex(double(img(:,:,2)), inv_indxs);
      img_g = reshape(img_g, size(X));
      img_b = bilinear_mex(double(img(:,:,3)), inv_indxs);
      img_b = reshape(img_b, size(X));

      img = cat(3, img_r, img_g, img_b);
    else
      img = bilinear_mex(double(img), inv_indxs);
      img = reshape(img, size(X));
    end
  end


  return;
end
