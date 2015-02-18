function [x, y, u, v, SnR] = matpiv_nfft(im1, im2, wins, overlap, thresh, mask)
% MATPIV_NFFT multiple passes with subpixel NFFT resolution from MatPIV.
%
%   [X, Y, U, V] = MATPIV_NFFT(IM1, IM2, WINSIZE, OVERLAP) performs PIV in
%   multiple passes with an interrogation window of size WINSIZE between IM1
%   and IM2 to eliminate the displacement bias, using continuous window shifting
%   in the Fourier domain [1]. OVERLAP defines the overlap between consecutive
%   windows.
%
%   [...] = MATPIV_NFFT(..., THRESH) defines the threshold THESH used to filter
%   out spurious vectors (both global and local filtering, default values: 3).
%
%   [...] = MATPIV_NFFT(..., MASK) provides a boolean MASK, the same size as IM1
%   and IM2, defining which positions of the images should be analyzed.
%
%   [..., SNR] = MATPIV_NFFT(...) returns in additon the SnR.
%
% References:
% [1] Liao, Q., & Cowen, E. a. "An efficient anti-aliasing spectral continuous 
% window shifting technique for PIV", Experiments in Fluids 38 (2005) 197–208.
%
% Quite heavly adapted by Simon Blanchoud (concatenated files, bug fix and speedup).
% 10.02.2015
%
% Copyright 1998-2011, Kristian Sveen, jks@math.uio.no/j.k.sveen@gmail.com
% for use with MatPIV 1.7
%
% Distributed under the terms of the Gnu General Public License manager
% Time stamp: 09:33, Mar 4 2011

  x = [];
  y = [];
  u = [];
  v = [];
  SnR = [];

  if (nargin < 4)
    disp('MatPIV requires at least 4 inputs. Aborting !')
    return;
  elseif (nargin < 5)
    thresh = 3;
    mask = [];
  elseif (nargin < 6)
    mask = [];
  end

  if (islogical(thresh))
    tmp = mask;
    mask = thresh;
    thresh = tmp;
  end

  im1 = double(im1);
  im2 = double(im2);
  imgsize = size(im1);

  if (isempty(mask))
    mask = true(imgsize);
  end

  if (isempty(thresh))
    thresh = 3;
  end

  snr_thresh = 1;

  if (any(imgsize ~= size(im2) | imgsize ~= size(mask)))
    disp('Images must have consistent sizes. Aborting !')
    return;
  end

  if size(wins,1)==1
      if size(wins,2)==1
          wins=[wins, wins];
      end
      wins=[wins; wins/2];
  end
  wins = ceil(wins);
  wins = wins + rem(wins, 2);

  iter=size(wins,1);

  if (numel(thresh) ~= iter)
    thresh = [thresh(:); ones(iter-numel(thresh), 1) * thresh(end)];
  end

  datax = [];
  datay = [];

  for i=1:iter-1

      [x,y,datax,datay,win_mask] = remesh(imgsize, wins(i,:), overlap, x, y, datax, datay, mask);

      [datax,datay,snr]=firstpass(im1,im2,wins(i,:),x,y,datax,datay,win_mask);

      [datax,datay]=snrfilt(datax,datay,snr,snr_thresh);
      [datax,datay]=globfilt(x,y,datax,datay,thresh(i));
      [datax,datay]=localfilt(x,y,datax,datay,thresh(i),win_mask);

      if (all(size(datax) > 1))
        [datax,datay]=naninterp2(datax,datay,win_mask,x,y);
      end

      datax=round(datax);
      datay=round(datay);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Final pass. Gives displacement to subpixel accuracy.
  [x,y,datax,datay,win_mask] = remesh(imgsize, wins(end,:), overlap, x, y, datax, datay, mask);

  [datax,datay,SnR]=finalpass(im1,im2,wins(end,:),x,y,datax,datay,win_mask);

  [datax,datay]=snrfilt(datax,datay,snr,snr_thresh);
  [datax,datay]=globfilt(x,y,datax,datay,thresh(end));
  [datax,datay]=localfilt(x,y,datax,datay,thresh(end),win_mask);

  if (all(size(datax) > 1))
    [u,v]=naninterp2(datax,datay,win_mask,x,y);
  else
    u = datax;
    v = datay;
  end

  return;
end

function [xx,yy,datax,datay,win_maske] = remesh(imgsize, winsize, ol, prevx, prevy, datax, datay, maske)

  M=winsize(1);
  N=winsize(2);
  x=[1:((1-ol)*M):imgsize(2)-M+1];
  y=[1:((1-ol)*N):imgsize(1)-N+1];

  nx=length(x);
  ny=length(y);

  xx=repmat(x+M/2,ny,1);
  yy=repmat((y+N/2).',1,nx);

  if (isempty(datax) || isempty(datay))
    datax = zeros(ny,nx);
    datay = zeros(ny,nx);
  elseif (numel(prevx) == 1)
    datax = ones(ny,nx) * datax;
    datay = ones(ny,nx) * datay;
  else
    datax = round(interp2(prevx,prevy,datax,xx,yy));
    datay = round(interp2(prevx,prevy,datay,xx,yy));
  end
  %datax(isnan(datax)) = 0;
  %datay(isnan(datay)) = 0;

  win_maske=(interp2(double(maske),xx,yy)>=0.5);

  return;
end

function [datax,datay,snr]=firstpass(A,B,N,xx,yy,idx,idy,maske)

% function [x,y,datax,datay]=firstpass(A,B,M,ol,idx,idy,maske)
%
% This function is used in conjunction with the MULTIPASS.M run-file.
% Inputs are allocated from within MULTIPASS.

% 1999 - 2005, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.5, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 14.41, 24 Mar 2005

  M=N(1); N=N(2); 
  [sy,sx]=size(A);

  x=round(xx(1,:)-M/2);
  y=round(yy(:,1)-N/2);

  nx=length(x);
  ny=length(y);

  datax=NaN(ny,nx);
  datay=NaN(ny,nx);
  snr=NaN(ny,nx);

  idx(isnan(idx)) = 0;
  idy(isnan(idy)) = 0;

  nelems=M*N;
  inelems = 1/nelems;
  if (nelems>1)
    nestim=nelems-1;
  else
    nestim=nelems;
  end
  inestim = 1/nestim;

  dev=20;
  tmpw = (1-cos(pi*(0:N-1)/(N-1)).^dev);tmpw2 = (1-cos(pi*(0:M-1)/(M-1)).^dev);
  W = tmpw'*tmpw2;

  mf = 2^nextpow2(M+N);

  full_sizes = [mf-1, 1/(mf-1)];
  sub_sizes = [N-4, 1/(N-4)];
  offset = [N M]/2 + 1;
  no_off = [0 0];

  for cj=1:ny
      jj = y(cj);
      for ci=1:nx
        ii = x(ci);
          if maske(cj,ci)
              if jj+idy(cj,ci)<1
                  idy(cj,ci)=1-jj;
              elseif jj+idy(cj,ci)>sy-N+1
                  idy(cj,ci)=sy-N+1-jj;
              end
              if ii+idx(cj,ci)<1
                  idx(cj,ci)=1-ii;
              elseif ii+idx(cj,ci)>sx-M+1
                  idx(cj,ci)=sx-M+1-ii;
              end

              E=A(jj:jj+N-1,ii:ii+M-1);   
              D2=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));

              mD2 = sum(D2(:)) * inelems;
              mE = sum(E(:)) * inelems;

              D2 = D2 - mD2;
              E = E - mE;

              stad1= sqrt(sum(E(:).^2) * inestim);
              stad2= sqrt(sum(D2(:).^2) * inestim);

              % use weights
              E = E.*W;
              F = D2.*W;

              E = E - sum(E(:)) * inelems;
              F = F - sum(F(:)) * inelems;

              % take zero-padded Fourier Transform
              at = fftn(E,[mf mf]);
              bt = fftn(conj(F(end:-1:1,end:-1:1)),[mf mf]);

              %%%%%%%%%%%%%%%%%%%%%% Calculate the normalized correlation:
              R = real(ifftn(bt.*at, 'nonsymmetric'));
              R=R(1:end-1,1:end-1);
              R=real(R)./(nelems*stad1*stad2);

              %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
              if full_sizes(1)==(N-1) || N < 5 || M < 5
                  [max_y1,max_x1,max_val]=getmax(R, full_sizes, no_off);
              else
                  [max_y1,max_x1,max_val]=getmax(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3), sub_sizes, offset);
              end

              datax(cj,ci)=(M-max_x1)+idx(cj,ci);
              datay(cj,ci)=(N-max_y1)+idy(cj,ci);

              snr(cj,ci) = mf * max_val.^2 / sum(R(:).^2);
          end
      end
  end

  return;
end

function [up,vp,SnR]=finalpass(A,B,N,xx,yy,idx,idy,maske)
% function [x,y,u,v,SnR,PeakHeight,brc]=finalpass_new(A,B,N,ol,idx,idy,Dt,mask)
%
% Provides the final pass to get the displacements with
% subpixel resolution. Uses sub-pixel displacements in Fourier Space
% following the paper by Qian and Cowen (2005, Experiments in Fluids).
%
%

% 1999 - 2011, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.7a, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 09:40, 4 Mar 2011

  M=N(1); N=N(2);
  [sy,sx]=size(A);

  x=round(xx(1,:)-M/2);
  y=round(yy(:,1)-N/2);

  nx=length(x);
  ny=length(y);

  up=NaN(ny,nx);
  vp=NaN(ny,nx);

  SnR=NaN(ny,nx);

  idx(isnan(idx)) = 0;
  idy(isnan(idy)) = 0;

  % Variables used for subpixel displacement in Fourier Domain
  min_res=0.005; % minimum residual to reach before breaking out of
  % sub-pixel iterations
  breakoutcounter=1; % count the iterations
  max_iterations=10; % max iterations before breaking out
  dev=20;
  tmpw = (1-cos(pi*(0:N-1)/(N-1)).^dev);tmpw2 = (1-cos(pi*(0:M-1)/(M-1)).^dev);
  W = tmpw'*tmpw2;

  %%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
  nelems=M*N;
  inelems = 1/nelems;
  if (nelems>1)
    nestim=nelems-1;
  else
    nestim=nelems;
  end
  inestim = 1/nestim;

  mf = 2^nextpow2(M+N);

  % window shift
  window_shift = zeros(mf);
  I=1:2*M; I(I>M)=I(I>M)-2*M; I=repmat(I,2*N,1); % used in the sub-pixel
  window_shift(1:2*N,1:2*M) = I;
  I=window_shift;
  J=(1:2*N)'; J(J>N)=J(J>N)-2*N; J=repmat(J,1,2*M);
  window_shift(1:2*N,1:2*M) = J;
  J=window_shift;

  iI = 1/(size(I,2));
  iJ = 1/(size(J,1));

  full_sizes = [mf-1, 1/(mf-1)];
  sub_sizes = [N-4, 1/(N-4)];
  alt_sizes = [N, 1/N];
  offset = [N M]/2 + 1;
  alt_off = [N M]/2 - 1;
  no_off = [0 0];

  for cj=1:ny
      jj = y(cj);
      for ci=1:nx
        ii = x(ci);
          if maske(cj,ci)
              if jj+idy(cj,ci)<1
                  idy(cj,ci)=1-jj;
              elseif jj+idy(cj,ci)>sy-N+1
                  idy(cj,ci)=sy-N+1-jj;
              end
              if ii+idx(cj,ci)<1
                  idx(cj,ci)=1-ii;
              elseif ii+idx(cj,ci)>sx-M+1
                  idx(cj,ci)=sx-M+1-ii;
              end
              E=A(jj:jj+N-1,ii:ii+M-1);
              D2=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));

              mD2 = sum(D2(:)) * inelems;
              mE = sum(E(:)) * inelems;

              D2 = D2 - mD2;
              E = E - mE;

              stad1= sqrt(sum(E(:).^2) * inestim);
              stad2= sqrt(sum(D2(:).^2) * inestim);

              ok1 = (stad1>eps);
              ok2 = (stad2>eps);

              if (ok1 && ok2)
                  % use weights
                  E = E.*W;
                  F = D2.*W;

                  E = E - sum(E(:)) * inelems;
                  F = F - sum(F(:)) * inelems;

                  % take zero-padded Fourier Transform
                  at = fftn(E,[mf mf]);
                  bt = fftn(conj(F(end:-1:1,end:-1:1)),[mf mf]);

                  %%%%%%%%%%%%%%%%%%%%%% Calculate the normalized correlation:
                  R = real(ifftn(bt.*at, 'nonsymmetric'));
                  R=R(1:end-1,1:end-1);
                  R=real(R)./(nelems*stad1*stad2);

                  %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
                  if full_sizes(1)==(N-1) || N < 5 || M < 5
                      [max_y1,max_x1,max_val]=getmax(R, full_sizes, no_off);
                  else
                      [max_y1,max_x1,max_val]=getmax(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3), sub_sizes, offset);
                  end

                  % loop on integer basis to make sure we've converged to
                  % +-0.5 pixels before entering subpixel shifting
                  stx=idx(cj,ci); sty=idy(cj,ci);

                  while max_x1~=M && max_y1~=N && ...
                          breakoutcounter<max_iterations &&...
                          jj+sty>0 && ii+stx>0 && ii+M-1+stx<=sx && jj+N-1+sty<=sy

                      D2=B(jj+sty:jj+N-1+sty,...
                          ii+stx:ii+M-1+stx);

                      F=(D2-sum(D2(:))*inelems).*W; F=F-sum(F(:))*inelems;
                      bt = fftn(conj(F(end:-1:1,end:-1:1)),[mf mf]);
                      R = real(ifftn(bt.*at, 'nonsymmetric'));

                      R=R(1:end-1,1:end-1);
                      R=real(R)./(nelems*stad1*stad2);
                      [max_y1,max_x1,max_val]=getmax(R, full_sizes, no_off);

                      stx=stx + (M-max_x1);
                      sty=sty + (N-max_y1);

                      breakoutcounter=breakoutcounter+1;
                  end

                  if breakoutcounter~=max_iterations
                      %update these only IF convergence was met, that is, we
                      %used less than max_iterations
                      idx(cj,ci)=stx; idy(cj,ci)=sty; 
                      breakoutcounter=1; % only reset if converged
                  end

                  %Only enter next bit if the peak is not located at the
                  %edges of the correlation plane
                  if max_x1~=1 && max_y1~=1 && max_x1~=mf-1 && max_y1~=mf-1
                      % 3-point peak fit using centroid, gaussian (default)
                      % or parabolic fit
                      [x0, y0]=intpeak(max_x1,max_y1,R(max_y1,max_x1),...
                          R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
                          R(max_y1-1,max_x1),R(max_y1+1,max_x1),M,N);
                      X0=x0; Y0=y0;
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % here we do the subpixel shifts in Fourier space
                      while (abs(x0)>min_res || abs(y0)>min_res) &&...
                              breakoutcounter<max_iterations

                          bt2=(exp(2*1i*pi*( (I*(X0)*iI) + ...
                              (J*(Y0)*iJ)))).*bt;
                          % At this point we could do some simple
                          % FFT-filtering like Todd and Liao suggests in
                          % their paper
                          % for example only keeping a certain number of
                          % Fourier components
                          % 26/3-05: This makes the peak a lot wider/broader:
                          % W2=1-weight('cosn',64,20)
                          % R=ifft2( bt2.*at.*W2 );
                          %
                          R=ifftn( bt2.*at, 'nonsymmetric');
                          R=real(R(1:end-1,1:end-1));
                          R(R<=0) = 1e-6;
                          R=real(R)./(nelems*stad1*stad2);

                          [dy,dx,max_val]=getmax(R, full_sizes, no_off);

                          X0=X0+(M-dx); Y0=Y0+(N -dy);
                          if dx>1 && dx<mf-1 && dy>1 && dy<mf-1
                              %only gaussian fit here
                              x0= -(log(R(dy,dx-1))-log(R(dy,dx+1)))/...
                                  (2*log(R(dy,dx-1))-4*log(R(dy,dx))+2*log(R(dy,dx+1)));
                              y0= -(log(R(dy-1,dx))-log(R(dy+1,dx)))/...
                                  (2*log(R(dy-1,dx))-4*log(R(dy,dx))+2*log(R(dy+1,dx)));

                              X0=X0-x0; Y0=Y0-y0;
                              breakoutcounter=breakoutcounter+1;
                          else
                              X0=nan; Y0=nan;
                              breakoutcounter=max_iterations;
                          end
                      end
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                      % Find the signal to Noise ratio
                      SnR(cj,ci) = mf * max_val.^2 / sum(R(:).^2);

                      %%%%%%%%%%%%%%%%%%%%%% Store the displacements, SnR and Peak Height.
                      up(cj,ci)=(-X0+idx(cj,ci));
                      vp(cj,ci)=(-Y0+idy(cj,ci));
                  end
              end
          end
          breakoutcounter=1;
      end
  end

  return;
end

function [u,v]=snrfilt(u,v,snr,thresh)

  bads = ~(snr >= thresh);

  u(bads) = NaN;
  v(bads) = NaN;

  return;
end

function [u,v]=globfilt(x,y,u,v,thresh)
% This function is a so called global histogram operator. It
% features a few slightly different ways of giving the maximum and
% minimum velocities allowed in your vector fields.
%
% There are two basically different methods used in GLOBFILT, The first
% uses a graphical input to specify the acceptance interval of velocity
% vectors. These are plotted in the (u,v) plane and the user should
% specify 4 points, using the left mouse button, that together form a 4
% sided polygonal region of acceptance. Use 'manual' as input parameter
% for this option.  Alternatively one can make an acception interval
% based on the standard deviations (x and y) of the measurement
% ensemble. This can be done in three different ways, namely 1) by
% specifying a factor (a number), 2) by specifying 'loop' or 3) by
% specifying a vector with the upper and lower velocitylimits. In the
% former case GLOBFILT uses the mean of the velocities plus/minus the
% number times the standard deviation as the limits for the acceptance
% area. In the second case GLOBFILT loops and lets the user
% interactively set the factor. This option often performs well if the
% vector field is not heavily contaminated with outliers. The third
% option is used by specifying an input vector [Umin Umax Vmin Vmax]
% which defines the upper and lower limits for the velocities.
%
% Additionally you can include 'interp' as a final input to
% interpolate the outliers found by GLOBFILT.
%
% See also MATPIV, SNRFILT, LOCALFILT, MASK

% 1999 -2014 copyright J.K.Sveen jks@math.uio.no
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
%
% For use with MatPIV 1.7 and later versions
% Distributed under the Gnu General Public License

  [xo, sx] = mymean(u(:));
  [yo, sy] = mymean(v(:));

  distx = (u(:) - xo).^2;
  disty = (v(:) - yo).^2;
  valids = (distx <= (thresh*sx)^2) & (disty <= (thresh*sy)^2);

  u(~valids)=NaN; v(~valids)=NaN;

  return;
end

function [u, v] = localfilt(x, y, u, v, threshold, maske)
% [NewU,NewV]=localfilt(x,y,u,v,threshold,mask)
%
% This function is a filter that will remove vectors that deviate from
% the median or the mean of their surrounding neighbors by the factor
% THRESHOLD times the standard deviation of the neighbors. 
%
% MASK should be applied to save calculation time. LOCALFILT is relatively 
% slow on large matrices and exlcuding "non-interesting" regions with MASK
% can increase speed in some cases.
%
% 1999 - 2002 , jks@math.uio.no
% For use with MatPIV 1.6
%
% Copyright J.K.Sveen (jks@math.uio.no)
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
% Distributed under the Gnu General Public License
%
% Time: 10:41, Jan 17 2002

  border = 2;
  valids = true(2*border + 1);
  valids(border+1, border+1) = false;
  nans = NaN(1, 2);

  histou = blockproc(u, [1 1], @blockstats, 'BorderSize', [border border], 'PadMethod', NaN, 'TrimBorder', false);
  histov = blockproc(v, [1 1], @blockstats, 'BorderSize', [border border], 'PadMethod', NaN, 'TrimBorder', false);

  medu = histou(:,1:2:end);
  stdu = threshold*histou(:,2:2:end);

  medv = histov(:,1:2:end);
  stdv = threshold*histov(:,2:2:end);

  bads = (u > medu + stdu | u < medu - stdu | ...
          v > medv + stdv | v < medv - stdv);

  u(bads)=NaN;
  v(bads)=NaN;

  return;

  function vals = blockstats(blk)

    if (maske(blk.location(1), blk.location(2)))
      data = blk.data(valids);
      vals = [nanmedian(data), nanstd(data)];
    else
      vals = nans;
    end

    return;
  end
end

function [si,sj,val] = getmax(mat, sizes, offset)

  [val, ind] = max(mat(:));
  si = rem(ind-1,sizes(1)) + 1;
  sj = (ind-si)*sizes(2) + 1 + offset(2);
  si = si + offset(1);

  return;
end

function [u,v]=naninterp2(u,v,mask,xx,yy)
% function [u,v]=naninterp2(u,v,mask,x,y)
%
% Interpolates NaN's in a vectorfield. Used by GLOBFILT,
% MEDIANFILT and VALIDATE. Sorts all spurious vectors based on the
% number of spurous neighbors to a point. Interpolation starts with
% the ones that have the least number of outliers in their
% neighborhood and loops until no NaN's are present in the field.
%


% 1999 - 2001, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.6, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 15:46, Oct 23, 2001

  % Replaced that function with the solution from John D'Errico
  % that is around 10x faster and works much more precisely,
  % in particular for smooth surfaces
  u = inpaint_nans(u, 2);
  u(~mask) = NaN;

  v = inpaint_nans(v, 2);
  v(~mask) = NaN;

  return;
end

function [x0,y0]=intpeak(x1,y1,R,Rxm1,Rxp1,Rym1,Ryp1,M,N)
% INTPEAK - interpolate correlation peaks in PIV
%
% function [x0,y0]=intpeak(x1,x2,x3,y1,y2,y3,method,N)
% METHOD = 
% 1 for centroid fit, 
% 2 for gaussian fit, 
% 3 for parabolic fit
% x1 and y1 are maximal values in respective directions.
% N is interrogation window size. N is either 1x1 or 1x2
%
% This is a subfunction to MATPIV

% Time stamp: 12:32, Apr. 14, 2004.
%
% Copyright 1998-2004, J. Kristian Sveen, 
% jks@math.uio.no/jks36@damtp.cam.ac.uk
% Dept of Mathmatics, University of Oslo/ 
% DAMTP, Univ. of Cambridge, UK
%
% Distributed under the GNU general public license.
%
% For use with MatPIV 1.6 and subsequent versions

  if any(([R Rxm1 Rxp1 Rym1 Ryp1])==0)
      % to avoid Log of Zero warnings
      x01=(((x1-1)*Rxm1)+(x1*R)+((x1+1)*Rxp1)) / (Rxm1+ R+Rxp1);
      y01=(((y1-1)*Rym1)+(y1*R)+((y1+1)*Ryp1)) / (Rym1+ R+Ryp1);
      x0=x01-(M);
      y0=y01-(N);
  else
      x01=x1 + ( (log(Rxm1)-log(Rxp1))/( (2*log(Rxm1))-(4*log(R))+(2*log(Rxp1))) );
      y01=y1 + ( (log(Rym1)-log(Ryp1))/( (2*log(Rym1))-(4*log(R))+(2*log(Ryp1))) );  
      x0=x01-(M);
      y0=y01-(N);
  end

  x0=real(x0);
  y0=real(y0);

  return;
end
