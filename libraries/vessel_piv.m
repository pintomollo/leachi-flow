function [x, y, u, v, snr] = vessel_piv(im1, im2, wins, overlap, mapping, centers, thresh, snr_thresh)
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
% Added a gaussian smoothing on the correlation matrix and SnR computation to reduce
% the impact of noise
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
  snr = [];

  if (nargin < 6)
    disp('Vessel PIV requires at least 6 inputs. Aborting !')
    return;
  elseif (nargin < 7)
    thresh = 3;
    snr_thresh = 1;
  elseif (nargin < 8)
    snr_thresh = 1;
  end

  min_n = 4;

  im1 = double(im1);
  im2 = double(im2);
  imgsize = size(im1);

  if (any(imgsize ~= size(im2) | imgsize ~= size(mapping)))
    disp('Images must have consistent sizes. Aborting !')
    return;
  end

  if size(wins,2)==1
    wins=[wins, wins];
  end
  if size(wins,1)==1
    wins=[wins; wins/2];
  end
  wins = ceil(wins);
  wins = wins + rem(wins, 2);

  iter=size(wins,1);

  if (numel(thresh) ~= iter)
    thresh = [thresh(:); ones(iter-numel(thresh), 1) * thresh(end)];
  end

  [vessel_params] = init_vessel(mapping, centers);

  datax = [];
  datay = [];
  group = [];

  out_dist = bwdist(logical(mapping));

  accum = zeros(max(mapping(:)), 3);

  for i=1:iter-1
      [x,y,group,datax,datay] = vessel_remesh(imgsize, wins(i,:), overlap, mapping, vessel_params, x, y, group, datax, datay);

      [datax,datay,snr]=firstpass(im1,im2,wins(i,:),x,y,datax,datay);

      [datax,datay]=snrfilt(datax,datay,snr,snr_thresh);
      [datax,datay]=maskfilt(x,y,datax,datay,out_dist,wins(i,:));
      [datax,datay]=globfilt(datax,datay,thresh(i));
      [datax,datay]=localfilt(datax,datay,thresh(i),group,min_n);

      [datax,datay]=groupinterp(datax,datay,group);

      [datax, datay, accum]=converge(datax,datay,group,accum);

      %{
      if (all(size(datax) > 1))
        [datax,datay]=naninterp2(datax,datay,win_mask,x,y);
      end

      goods = (~isnan(datax));
      if (any(goods(:)))
        accum(:,:,1) = accum(:,:,1) + datax;
        accum(:,:,2) = accum(:,:,2) + datay;
        accum(:,:,3) = accum(:,:,3) + goods;
      end
      %}

      datax=round(datax);
      datay=round(datay);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Final pass. Gives displacement to subpixel accuracy.
  [x,y,group,datax,datay] = vessel_remesh(imgsize, wins(end,:), overlap, mapping, vessel_params, x, y, group, datax, datay);

  %datax = round(accum(:,:,1) ./ accum(:,:,3));
  %datay = round(accum(:,:,2) ./ accum(:,:,3));

  [datax,datay,snr]=finalpass(im1,im2,wins(end,:),x,y,datax,datay);

  [datax,datay]=snrfilt(datax,datay,snr,snr_thresh);
  [datax,datay]=maskfilt(x,y,datax,datay,out_dist,wins(end,:));
  [datax,datay]=globfilt(datax,datay,thresh(end));
  [datax,datay]=localfilt(datax,datay,thresh(end),group,min_n);

  %{
  goods = ~isnan(datax);

  if (all(size(datax) > 1))
    [u,v]=naninterp2(datax,datay,win_mask,x,y);
  else
    u = datax;
    v = datay;
  end

  [u, v]=converge(u,v,accum);
  u(goods) = datax(goods);
  v(goods) = datay(goods);
  %}
  [datax, datay]=groupinterp(datax,datay,group);

  [u, v]=converge(datax,datay,group,accum);

  dist = ((u - datax).^2 + (v - datay).^2);
  refined = (dist < 1);

  u(refined) = datax(refined);
  v(refined) = datay(refined);

  return;
end

function vessel_params = init_vessel(mapping, centers)

  x1 = centers(1:3:end, 1);
  x2 = centers(2:3:end, 1);
  y1 = centers(1:3:end, 2);
  y2 = centers(2:3:end, 2);

  vects = [x2-x1 y2-y1];
  middle = ([x1 y1] + [x2 y2])/2;
  lens = sqrt(sum(vects.^2, 2));
  vects = bsxfun(@rdivide, vects, lens);
  widths = zeros(size(lens));
  angl = zeros(size(lens));

  nbranches = length(x1);
  max_dist = sqrt(sum(size(mapping).^2));
  y = [0:max_dist];
  center_indx = length(y);
  y = [-y(end:-1:2) y];

  for i=1:nbranches
    curr_branch = double(mapping == i);
    x = [0:lens(i)];

    [X,Y] = meshgrid(x,y);
    angl(i) = -atan2(vects(i,2), vects(i,1));
    rot = [cos(angl(i)) -sin(angl(i)); sin(angl(i)) cos(angl(i))];

    pts = [X(:) Y(:)] * rot;
    pts = bsxfun(@plus, pts, [x1(i) y1(i)]);

    proj = bilinear_mex(curr_branch, pts);
    proj = reshape(proj, size(X));

    proj = proj(center_indx:end,:) + proj(center_indx:-1:1,:);
    rad = nansum(proj,2);
    widths(i) = find(rad==0, 1, 'first');
  end

  vessel_params = [middle angl lens/2 widths];

  return;
end

function [xx,yy,gg,datax,datay] = vessel_remesh(imgsize, winsize, ol, mapping, vessel, prevx, prevy, prevg, datax, datay)

  M=winsize(1);
  N=winsize(2);

  rim = ceil(winsize/2);

  xx=[];
  yy=[];
  gg=[];

  origs = vessel(:,1:2);
  angls = vessel(:,3);
  lens = vessel(:,4);
  widths = vessel(:,5);

  nbranches = size(vessel, 1);
  for i=1:nbranches
    x=[0:((1-ol)*M):lens(i)];
    y=[0:((1-ol)*N):widths(i)];

    x = [-x(end:-1:2) x];
    y = [-y(end:-1:2) y];

    [X,Y] = meshgrid(x,y);
    rot = [cos(angls(i)) -sin(angls(i)); sin(angls(i)) cos(angls(i))];

    pts = [X(:) Y(:)] * rot;
    pts = round(bsxfun(@plus, pts, origs(i,:)));

    valids = all(bsxfun(@gt, pts, rim) & bsxfun(@le, pts, imgsize([2 1])-rim), 2);
    pts = pts(valids, :);

    groups = mapping((pts(:,1)-1)*imgsize(1) + pts(:,2));
    goods = (groups == i);

    xx = [xx; pts(goods, 1)];
    yy = [yy; pts(goods, 2)];
    gg = [gg; groups(goods)];
  end

  new_size = [length(xx) 1];

  valids = (~isnan(datax) & ~isnan(datay));

  if (isempty(datax) || isempty(datay) || ~any(valids) || isempty(xx))
    datax = zeros(new_size);
    datay = zeros(new_size);
  %elseif (numel(prevx) == 1)
  %  datax = ones(new_size) * datax;
  %  datay = ones(new_size) * datay;
  else

    dists = bsxfun(@minus, xx, prevx(valids).').^2 + ...
            bsxfun(@minus, yy, prevy(valids).').^2;

    sames = bsxfun(@eq, gg, prevg(valids).');

    dists(~sames) = Inf;

    weights = exp(-bsxfun(@rdivide, dists, 2*widths(gg).^2));
    norms = sum(weights, 2);

    weights = bsxfun(@rdivide, weights, norms);

    datax = sum(bsxfun(@times, weights, datax(valids).'), 2);
    datay = sum(bsxfun(@times, weights, datay(valids).'), 2);

    datax = round(datax);
    datay = round(datay);
  end

  return;
end

function [datax,datay,snr]=firstpass(A,B,N,xx,yy,idx,idy)
%
% This function is used in conjunction with the MULTIPASS.M run-file.
% Inputs are allocated from within MULTIPASS.

% 1999 - 2005, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.5, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 14.41, 24 Mar 2005

  M=N(1); N=N(2); 
  [sy,sx]=size(A);

  x=round(xx-M/2);
  y=round(yy-N/2);

  n=length(x);

  datax=NaN(n,1);
  datay=NaN(n,1);
  snr=NaN(n,1);

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

  for i=1:n
    jj = y(i);
    ii = x(i);

    if jj+idy(i)<1
        idy(i)=1-jj;
    elseif jj+idy(i)>sy-N+1
        idy(i)=sy-N+1-jj;
    end
    if ii+idx(i)<1
        idx(i)=1-ii;
    elseif ii+idx(i)>sx-M+1
        idx(i)=sx-M+1-ii;
    end

    E=A(jj:jj+N-1,ii:ii+M-1);
    D2=B(jj+idy(i):jj+N-1+idy(i),ii+idx(i):ii+M-1+idx(i));

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
        R_peak=gaussian_mex(R, 2);

        %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
        if full_sizes(1)==(N-1) || N < 5 || M < 5
            [max_y1,max_x1,max_val]=getmax(R_peak, full_sizes, no_off);
        else
            [max_y1,max_x1,max_val]=getmax(R_peak(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3), sub_sizes, offset);
        end

        datax(i)=(M-max_x1)+idx(i);
        datay(i)=(N-max_y1)+idy(i);

        snr(i) = mf * max_val.^2 / sum(R_peak(:).^2);
    end
  end

  return;
end

function [up,vp,SnR]=finalpass(A,B,N,xx,yy,idx,idy)
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

  x=round(xx-M/2);
  y=round(yy-N/2);

  n=length(x);

  up=NaN(n,1);
  vp=NaN(n,1);

  SnR=NaN(n,1);

  idx(isnan(idx)) = 0;
  idy(isnan(idy)) = 0;

  % Variables used for subpixel displacement in Fourier Domain
  min_res=0.005; % minimum residual to reach before breaking out of
  % sub-pixel iterations
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
  window_shift=[1:mf/2 [mf/2+1:mf]-mf];
  I=repmat(window_shift, mf, 1);
  J=repmat(window_shift.', 1, mf);

  iI = 1/(size(I,2));
  iJ = 1/(size(J,1));

  full_sizes = [mf-1, 1/(mf-1)];
  sub_sizes = [N-4, 1/(N-4)];
  alt_sizes = [N, 1/N];
  offset = [N M]/2 + 1;
  alt_off = [N M]/2 - 1;
  no_off = [0 0];

  for i=1:n
    jj = y(i);
    ii = x(i);

    breakoutcounter=1; % count the iterations

    if jj+idy(i)<1
        idy(i)=1-jj;
    elseif jj+idy(i)>sy-N+1
        idy(i)=sy-N+1-jj;
    end
    if ii+idx(i)<1
        idx(i)=1-ii;
    elseif ii+idx(i)>sx-M+1
        idx(i)=sx-M+1-ii;
    end
    E=A(jj:jj+N-1,ii:ii+M-1);
    D2=B(jj+idy(i):jj+N-1+idy(i),ii+idx(i):ii+M-1+idx(i));

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
        R = real(ifftn(bt.*at, [mf mf], 'nonsymmetric'));
        R=R(1:end-1,1:end-1);
        R=real(R)./(nelems*stad1*stad2);
        R_peak=gaussian_mex(R, 2);

        %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
        if full_sizes(1)==(N-1) || N < 5 || M < 5
            [max_y1,max_x1,max_val]=getmax(R_peak, full_sizes, no_off);
        else
            [max_y1,max_x1,max_val]=getmax(R_peak(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3), sub_sizes, offset);
        end

        % loop on integer basis to make sure we've converged to
        % +-0.5 pixels before entering subpixel shifting
        stx=(M-max_x1)+idx(i);
        sty=(N-max_y1)+idy(i);

        while (max_x1~=M || max_y1~=N) && ...
                breakoutcounter<max_iterations &&...
                jj+sty>0 && ii+stx>0 && ii+M-1+stx<=sx && jj+N-1+sty<=sy

            D2=B(jj+sty:jj+N-1+sty,...
                ii+stx:ii+M-1+stx);

            F=(D2-sum(D2(:))*inelems).*W; F=F-sum(F(:))*inelems;
            bt = fftn(conj(F(end:-1:1,end:-1:1)),[mf mf]);
            R = real(ifftn(bt.*at, 'nonsymmetric'));

            R=R(1:end-1,1:end-1);
            R=real(R)./(nelems*stad1*stad2);
            R_peak=gaussian_mex(R, 2);
            [max_y1,max_x1,max_val]=getmax(R_peak, full_sizes, no_off);

            stx=stx + (M-max_x1);
            sty=sty + (N-max_y1);

            breakoutcounter=breakoutcounter+1;
        end

        if (breakoutcounter<max_iterations || (max_x1==M && max_y1==N))
            %update these only IF convergence was met, that is, we
            %used less than max_iterations
            idx(i)=stx; idy(i)=sty; 
            breakoutcounter=1; % only reset if converged
        else
            continue;
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

            dy = max_y1;
            dx = max_x1;

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

                if dx>1 && dx<mf-1 && dy>1 && dy<mf-1
                    %only gaussian fit here
                    x0=-(log(R(dy,dx-1))-log(R(dy,dx+1)))/...
                        (2*log(R(dy,dx-1))-4*log(R(dy,dx))+2*log(R(dy,dx+1)));
                    y0=-(log(R(dy-1,dx))-log(R(dy+1,dx)))/...
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
            R_peak=gaussian_mex(R, 2);
            SnR(i) = mf * R_peak(dy,dx).^2 / sum(R_peak(:).^2);

            %%%%%%%%%%%%%%%%%%%%%% Store the displacements, SnR and Peak Height.
            up(i)=(-X0+idx(i));
            vp(i)=(-Y0+idy(i));
        end
    end
  end

  return;
end

function [u,v]=snrfilt(u,v,snr,thresh)
% Removes data that have a poor SnR

  bads = (snr < thresh);

  u(bads) = NaN;
  v(bads) = NaN;

  return;
end

function [u,v]=maskfilt(x,y,u,v,mask,win_size)

  bads = isnan(u(:));
  u(bads) = 0;
  v(bads) = 0;

  valids = bilinear_mex(double(mask), x(:)+u(:), y(:)+v(:));
  valids = (valids < 0.5*mean(win_size));

  u(bads | ~valids) = NaN;
  v(bads | ~valids) = NaN;

  return;
end

function [u,v]=globfilt(u,v,thresh)
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

function [u, v] = localfilt(u, v, threshold, groups, min_n)
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

  [means, stds] = mymean([u v], 1, groups);

  [junk, indxs, reverse] = unique(groups);
  ngroups = diff([indxs; length(groups)+1]);

  means(ngroups < min_n, :) = NaN;

  means = means(reverse,:);
  stds = threshold*stds(reverse,:);

  bads = (u > means(:,1) + stds(:,1) | u < means(:,1) - stds(:,1) | ...
          v > means(:,2) + stds(:,2) | v < means(:,2) - stds(:,2));

  u(bads)=NaN;
  v(bads)=NaN;

  return;
end

function [datax,datay,full_accum] = converge(datax, datay, groups, accum)

  if (nargout > 2)
    full_accum = accum;

    [means] = mymean([datax datay], 1, groups);

    [curr_groups, indxs] = unique(groups);
    ngroups = diff([indxs; length(groups)+1]);

    means = bsxfun(@times, means, ngroups);

    goods = ~any(isnan(means), 2);
    if (any(goods))
      full_accum(curr_groups(goods), :) = full_accum(curr_groups(goods), :) + [means(goods,:) ngroups(goods)];
    end
  end

  accum = accum(groups, :);

  empties = isnan(datax);

  tmp = datax;
  tmp(empties) = 0;

  datax = (tmp + accum(:,1)) ./ (accum(:,3)+(~isnan(datax)));

  tmp = datay;
  tmp(empties) = 0;

  datay = (tmp + accum(:,2)) ./ (accum(:,3)+(~isnan(datay)));

  return;
end

function [si,sj,val] = getmax(mat, sizes, offset)

  [val, ind] = max(mat(:));
  si = rem(ind-1,sizes(1)) + 1;
  sj = (ind-si)*sizes(2) + 1 + offset(2);
  si = si + offset(1);

  return;
end

function [u,v]=groupinterp(u,v,groups)
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

  [means] = mymean([u v], 1, groups);

  [junk, indxs, reverse] = unique(groups);

  means = means(reverse,:);

  bads = (isnan(u) | isnan(v));

  u(bads) = means(bads, 1);
  v(bads) = means(bads, 2);

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
