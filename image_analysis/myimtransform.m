function [img] = myimtransform(img, modello, GLOBAL, final_size, origin)
% Just calling the MEX function that will do the actual warping: imwarp_mex

  flag_proj = (modello(1)~='a');
  proj_mat = inv(GLOBAL');
  type = class(img);
  img = imwarp_mex(double(img), proj_mat, flag_proj, final_size, origin);
  img = cast(img, type);

  return;
end
