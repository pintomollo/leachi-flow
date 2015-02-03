function [x,y,u,v,SnR]=matpiv_nfft(im1,im2,wins,overlap,sensit,maske,iter)
% MULTIPASSX - multiple passes 
% function [x,y,u,v,snr,pkh]=multipassx(im1,im2,winsize,time,...
% overlap,sensit,maske,numofiterations,ustart,vstart)
%
% PIV in multiple passes to eliminate the displacement bias.
% Utilizes the increase in S/N by  halving the size of the 
% interrogation windows after the first pass.
% Sub-function to MATPIV.
%
% This function uses continuous window shifting in the Fourier domain as
% published by Qian and Cowen (2005).
%
% See also: 
%          MATPIV, SNRFILT, LOCALFILT, GLOBFILT, DEFINEWOCO
  %
% Copyright 1998-2011, Kristian Sveen, jks@math.uio.no/j.k.sveen@gmail.com 
% for use with MatPIV 1.7
% Distributed under the terms of the Gnu General Public License manager
% Time stamp: 09:33, Mar 4 2011

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image read
  A=double(im1); B=double(im2);

  %%%%%%%%% First pass to estimate displacement in integer values:
  if nargin==4
      sensit=3;maske=true(size(im1)); iter=3;
  end
  if (isempty(maske))
    maske=true(size(im1));
  end
  [sy,sx]=size(A);

  if size(wins,1)==1
      if size(wins,2)==1
          wins=[wins, wins];
      end
      for jj=1:iter-2
          wins=[wins; wins(end,:)/2];
      end
  end
  iter=size(wins,1);

  for i=1:iter-1
      disp(['* Pass No: ',num2str(i)])
      if i==1
          [x,y,datax,datay]=firstpass(A,B,wins(i,:),overlap,[],[],maske);
      else
          [x,y,datax,datay]=firstpass(A,B,wins(i,:),overlap,datax,datay,maske);
      end
      win_maske=(interp2(double(maske),x,y)>=0.5);

      [datax,datay]=globfilt(x,y,datax,datay,3);
      [datax,datay]=localfilt(x,y,datax,datay,sensit,3,win_maske);
      [datax,datay]=naninterp2(datax,datay,win_maske,x,y);
      datax=floor(datax); datay=floor(datay);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % expand the velocity data to twice the original size
      %if i~=iter-1 
      %if wins(i,1)~=wins(i+1,1)    
      %  X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
      %  XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
      %else
      %  XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2; X=XI;
      %end
      %if wins(i,2)~=wins(i+1,2)
      %  Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
      %  YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2; 
      %else
      %  YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2; Y=YI;
      %end
      steps=wins(i+1,:) ./ wins(i,:);
      XI=[steps(1):steps(1):(size(datax, 2) + 1 - steps(1))];
      YI=[steps(2):steps(2):(size(datay, 1) + 1 - steps(2))];

      disp('   Expanding velocity-field for next pass')
      %datax=round(interp2(X,Y',datax,XI,YI'));
      %datay=round(interp2(X,Y',datay,XI,YI'));
      datax=round(interp2(datax,XI,YI'));
      datay=round(interp2(datay,XI,YI'));
      win_maske=(interp2(double(win_maske),XI,YI')>=0.5);
      [datax,datay]=naninterp2(datax,datay,win_maske,...
                              repmat(XI,size(datax,1),1),repmat(YI',1,size(datax,2)));

      datax=round(datax); datay=round(datay);
      %end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Final pass. Gives displacement to subpixel accuracy.
  disp('* Final Pass')
  [x,y,u,v,SnR]=finalpass_new(A,B,wins(end,:),overlap,round(datax),...
                              round(datay),maske);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(' ')

  return;
end

function [xx,yy,datax,datay]=firstpass(A,B,N,ol,idx,idy,maske)

% function [x,y,datax,datay]=firstpass(A,B,M,ol,idx,idy,maske)
%
% This function is used in conjunction with the MULTIPASS.M run-file.
% Inputs are allocated from within MULTIPASS.

% 1999 - 2005, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.5, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 14.41, 24 Mar 2005

  M=N(1); N=N(2); 
  overlap=ol; [sy,sx]=size(A);
  if isempty(idx) || isempty(idy)
      idx=zeros(floor(sy/(N*(1-ol))),floor(sx/(M*(1-ol))));
      idy=zeros(floor(sy/(N*(1-ol))),floor(sx/(M*(1-ol))));
  end
  xx=zeros(ceil((size(A,1)-N)/((1-overlap)*N))+1, ...
      ceil((size(A,2)-M)/((1-overlap)*M)) +1);
  yy=xx; datax=xx; datay=xx; 
  % change . october 2001, weight matrix added.
  % W=weight('cosn',[M N],100);
  %if nargin==7, 
  %    if ~isempty(maske)
  %        IN=zeros(size(maske(1).msk));
  %        for i=1:length(maske)
  %            IN=IN+double(maske(i).msk);
  %        end
  %    else 
  %        IN=zeros(size(A)); 
  %    end,
  %elseif nargin<8
  %    IN=zeros(size(A)); 
  %end

  padn = 2^nextpow2(2*N);
  padm = 2^nextpow2(2*M);

  cj=1;tic
  for jj=1:((1-ol)*N):sy-N+1
      ci=1;
      for ii=1:((1-ol)*M):sx-M+1 

          %if IN(jj+N/2,ii+M/2)~=1 
          if maske(jj+N/2,ii+M/2)

              if isnan(idx(cj,ci))
                  idx(cj,ci)=0;
              end
              if isnan(idy(cj,ci))
                  idy(cj,ci)=0;
              end
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

              C=A(jj:jj+N-1,ii:ii+M-1);   
              D=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));
  %             D=eval('B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci) )',...
  %                 'B(jj:jj+N-1,ii:ii+M-1)');
  %             if D==B(jj:jj+N-1,ii:ii+M-1)
  %                 idx(cj,ci)=0; idy(cj,ci)=0;
  %             end
              C=C-mean(C(:)); D=D-mean(D(:)); %C(C<0)=0; D(D<0)=0;
              stad1=std(C(:)); stad2=std(D(:)); 
              % Apply weight function by uncommenting below
              %C=C.*W; %D=D.*W;
              %
              if stad1==0, stad1=nan;end
              if stad2==0, stad2=nan; end

              %%%%%%%%%%%%%%%%%%%%%%%Calculate the normalized correlation:   
              R=xcorrf2(C,D,padn,padm)/(N*M*stad1*stad2);
              %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
              if size(R,1)==(N-1)
                [max_y1,max_x1]=find(R==max(R(:)));
              else
                [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3))));
              end

              if length(max_x1)>1
                max_x1=round(sum(max_x1.*(1:length(max_x1))')./sum(max_x1));
                max_y1=round(sum(max_y1.*(1:length(max_y1))')./sum(max_y1));
              elseif isempty(max_x1)
                idx(cj,ci)=nan; idy(cj,ci)=nan; max_x1=nan; max_y1=nan;
              end
              %%%%%%%%%%%%%%%%%%%%%% Store the displacements in variable datax/datay
              datax(cj,ci)=-(max_x1-(M))+idx(cj,ci);
              datay(cj,ci)=-(max_y1-(N))+idy(cj,ci);
              xx(cj,ci)=ii+M/2; yy(cj,ci)=jj+N/2;
              ci=ci+1;
          else
              xx(cj,ci)=ii+M/2; yy(cj,ci)=jj+N/2;
              datax(cj,ci)=NaN; datay(cj,ci)=NaN; ci=ci+1;
          end  
      end
      fprintf('\r No. of vectors: %d', ((cj-1)*(ci)+ci-1)-sum(isnan(datax(:))))
      fprintf(' , Seconds taken: %f', toc);
      cj=cj+1;
  end
  disp('.')

  return;
end

% now we inline the function XCORRF2 to shave off some time.
function c = xcorrf2(a,b,padn,padm)
%  c = xcorrf2(a,b)
%   Two-dimensional cross-correlation using Fourier transforms.
%       XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
%       XCORRF2(A) is the autocorrelation function.
%       This routine is functionally equivalent to xcorr2 but usually faster.
%       See also XCORR2.

%       Author(s): R. Johnson
%       $Revision: 1.0 $  $Date: 1995/11/27 $

  %if nargin==2
  %  pad='yes';
  %end

  [N,M] = size(a);

  %[ma,na] = size(a);
  %if nargin == 1
  %  %       for autocorrelation
  %  b = a;
  %end
  %[mb,nb] = size(b);
  %       make reverse conjugate of one array
  b = conj(b(end:-1:1,end:-1:1));
  %if strcmp(pad,'yes');  
    %       use power of 2 transform lengths
  %  mf = 2^nextpow2(ma+mb);
  %  nf = 2^nextpow2(na+nb);

    at = fftn(b,[padn,padm]);
    bt = fftn(a,[padn,padm]);
  %elseif strcmp(pad,'no');
  %  at = fft2(b);
  %  bt = fft2(a);
  %else
  %  disp('Wrong input to XCORRF2'); return
  %end

  %       multiply transforms then inverse transform
  c = ifftn(at.*bt, 'nonsymmetric');
  %       make real output for real input

  %if ~any(any(imag(a))) && ~any(any(imag(b)))
    c = real(c);
  %end
  %if strcmp(pad,'yes');
    %  trim to standard size
    %c(ma+mb:mf,:) = [];
    %c(:,na+nb:nf) = [];
    c=c(1:2*N,1:2*M);
  %  elseif strcmp(pad,'no');
  %  c=(c(1:end-1,1:end-1));

%    c(ma+mb:mf,:) = [];
%    c(:,na+nb:nf) = [];
  %end

  return;
end

function [xp,yp,up,vp,SnR,Pkh,brc]=finalpass_new(A,B,N,ol,idx,idy,maske)
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

  %  if length(N)==1
  %    M=N;
  %else
      M=N(1); N=N(2);
  %end
  cj=1;
  [sy,sx]=size(A);

  % Allocate space for matrixes
  xp=zeros(ceil((size(A,1)-N)/((1-ol)*N))+1, ...
      ceil((size(A,2)-M)/((1-ol)*M))+1);
  yp=xp; up=xp; vp=xp; brc=xp; SnR=xp; Pkh=xp;

  % Variables used for subpixel displacement in Fourier Domain
  min_res=0.005; % minimum residual to reach before breaking out of
  % sub-pixel iterations
  breakoutcounter=1; % count the iterations
  max_iterations=10; % max iterations before breaking out
  I=1:2*M; I(I>M)=I(I>M)-2*M; I=repmat(I,2*N,1); % used in the sub-pixel
  % window shift
  J=(1:2*N)'; J(J>N)=J(J>N)-2*N; J=repmat(J,1,2*M);
  %W=weight('cosn',[M,N],20); % weights used in the sub-pixel window shift
  dev=20;
  tmpw = (1-cos(pi*(0:M-1)/(M-1)).^dev);tmpw2 = (1-cos(pi*(0:N-1)/(N-1)).^dev);
  W = tmpw'*tmpw2;
  %W2=1-weight('cosn',2*[M N],20); %weight for FFT filtering

  %if nargin==7, 
  %    if ~isempty(maske)
  %        IN=zeros(size(maske(1).msk));
  %        for ii=1:length(maske)
  %            IN=IN+double(maske(ii).msk);
  %        end
  %    else IN=zeros(size(A)); 
  %    end, 
  %end

  fprintf([' Continuous windows shifting in Fourier space\n', ...
      '  Local iterations applied\n',...
      '  - Using ',num2str(M),'*',num2str(N),...
      ' interrogation windows! \n'])
  %%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
  nelems=M*N;
  if (nelems>1)
    nestim=nelems-1;
  else
    nestim=nelems;
  end

  mf = 2^nextpow2(M+N);
  nf = mf;

  tic
  for jj=1:((1-ol)*N):sy-N+1
      ci=1;
      for ii=1:((1-ol)*M):sx-M+1
          %if IN(jj+N/2,ii+M/2)~=1
          if maske(jj+N/2,ii+M/2)
              if isnan(idx(cj,ci))
                  idx(cj,ci)=0;
              end
              if isnan(idy(cj,ci))
                  idy(cj,ci)=0;
              end
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
              D2=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));
              E=A(jj:jj+N-1,ii:ii+M-1);

              mD2 = sum(D2(:)) / nelems;
              mE = sum(E(:)) / nelems;
              
              D2 = D2 - mD2;
              E = E - mE;

              stad1= sqrt(sum(E(:).^2) / nestim);
              stad2= sqrt(sum(D2(:).^2) / nestim);

              %stad1=std(E(:));
              %stad2=std(D2(:));
              if stad1==0, stad1=1; end
              if stad2==0, stad2=1; end

              % use weights
              E = E.*W;
              F = D2.*W;

              E = E - sum(E(:)) / nelems;
              F = F - sum(F(:)) / nelems;
              %E=(E-mean(E(:))).*W;
              %F=(D2-mean(D2(:))).*W;
              %E=E-mean(E(:));
              %F=F-mean(F(:));

              % take zero-padded Fourier Transform
              %mf = 2^nextpow2(M+N);
              %nf = mf;
              at = fftn(E,[nf mf]);
              bt = fftn(conj(F(end:-1:1,end:-1:1)),[nf mf]);

              %at = fft2(E,nf,mf);
              %bt = fft2(conj(F(end:-1:1,end:-1:1)),nf,mf);
              %no zero-padding version - need to change I and J if this
              %is uncommented
              %at = fft2(E);
              %bt = fft2(F);
              %%%%%%%%%%%%%%%%%%%%%% Calculate the normalized correlation:
              R = real(ifftn(bt.*at, 'nonsymmetric'));
              %R=real(ifft2(bt.*at));
              R(end,:)=[]; R(:,end)=[];
              R=real(R)./(nelems*stad1*stad2);
              %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
              %%%%%%%%%%%%%%%%%%%%%% _IF_ the standard deviation is NOT NaN.
              if all(~isnan(R(:))) && ~all(R(:)==0)  %~isnan(stad1) & ~isnan(stad2)
                  if size(R,1)==(N-1)
                      [max_y1,max_x1]=find(R==max(R(:)));

                  else
                      [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,...
                          0.5*M+2:1.5*M-3))));
                  end
                  if length(max_x1)>1
                      max_x1=round(sum(max_x1.^2)./sum(max_x1));
                      max_y1=round(sum(max_y1.^2)./sum(max_y1));
                  end

                  % loop on integer basis to make sure we've converged to
                  % +-0.5 pixels before entering subpixel shifting
                  stx=idx(cj,ci); sty=idy(cj,ci);

                  while max_x1~=M && max_y1~=N && ...
                          breakoutcounter<max_iterations &&...
                          jj+sty>0 && ii+stx>0 && ii+M-1+stx<=sx && jj+N-1+sty<=sy
                      D2=B(jj+sty:jj+N-1+sty,...
                          ii+stx:ii+M-1+stx);

                      F=(D2-sum(D2(:))/nelems).*W; F=F-sum(F(:))/nelems;
                      bt = fftn(conj(F(end:-1:1,end:-1:1)),[nf mf]);
                      R = real(ifftn(bt.*at, 'nonsymmetric'));

                      %F=(D2-mean(D2(:))).*W; F=F-mean(F(:));
                      %bt = fft2(conj(F(end:-1:1,end:-1:1)),nf,mf);
                      %R=ifft2(bt.*at);
                      R(end,:)=[]; R(:,end)=[];
                      R=real(R)./(nelems*stad1*stad2);
                      [max_y1,max_x1]=find(R==max(R(:)));
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
                  if max_x1~=1 && max_y1~=1 && max_x1~=M-1 && max_y1~=N-1
                      % 3-point peak fit using centroid, gaussian (default)
                      % or parabolic fit
                      [x0 y0]=intpeak(max_x1,max_y1,R(max_y1,max_x1),...
                          R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
                          R(max_y1-1,max_x1),R(max_y1+1,max_x1),M,N);
                      X0=x0; Y0=y0;
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % here we do the subpixel shifts in Fourier space
                      while (abs(x0)>min_res || abs(y0)>min_res) &&...
                              breakoutcounter<max_iterations
                          bt2=(exp(2*1i*pi*( (I*(X0)/size(I,2)) + ...
                              (J*(Y0)/size(J,1))))).*bt;
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
                          %R=ifft2( bt2.*at );
                          R(end,:)=[]; R(:,end)=[];
                          R=real(R)./(N*M*stad1*stad2);

                          [dy,dx]=find(R==max(R(:)));
                          X0=X0+(M-dx); Y0=Y0+(N -dy);
                          if dx>1 && dx<2*M-1 && dy>1 && dy<2*N-1
                              %only gaussian fit here
                              x0= -(log(R(dy,dx-1))-log(R(dy,dx+1)))/...
                                  (2*log(R(dy,dx-1))-4*log(R(dy,dx))+2*log(R(dy,dx+1)));
                              y0= -(log(R(dy-1,dx))-log(R(dy+1,dx)))/...
                                  (2*log(R(dy-1,dx))-4*log(R(dy,dx))+2*log(R(dy+1,dx)));

                              X0=X0-x0; Y0=Y0-y0;
                              breakoutcounter=breakoutcounter+1;
                              %if ~isreal(x0) | ~isreal(y0)
                              %  disp([num2str([ci,cj]),', ',...
                              %	    num2str([stad1 stad2 std(A(:))])])
                              %    end
                          else
                              X0=nan; Y0=nan;
                              breakoutcounter=16;
                          end
                      end
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                      % Find the signal to Noise ratio
                      R2=R;
                      try
                          R2(max_y1-3:max_y1+3,max_x1-3:max_x1+3)=NaN;
                      catch
                          R2(max_y1-1:max_y1+1,max_x1-1:max_x1+1)=NaN;
                      end
                      if size(R,1)==(N-1)
                          [p2_y2,p2_x2]=find(R2==max(R2(:)));                        
                      else
                          [p2_y2,p2_x2]=find(R2==max(max(R2(0.5*N:1.5*N-1,0.5*M:1.5*M-1))));
                      end
                      if length(p2_x2)>1
                          p2_x2=p2_x2(round(length(p2_x2)/2));
                          p2_y2=p2_y2(round(length(p2_y2)/2));
                      elseif isempty(p2_x2)

                      end
                      % signal to noise:
                      snr=R(max_y1,max_x1)/R2(p2_y2,p2_x2);
                      % signal to mean:
                      %snr=R(max_y1,max_x1)/mean(R(:));
                      % signal to median:
                      %snr=R(max_y1,max_x1)/median(median(R(0.5*N+2:1.5*N-3,...
                      %    0.5*M+2:1.5*M-3)));

                      %%%%%%%%%%%%%%%%%%%%%% Store the displacements, SnR and Peak Height.
                      up(cj,ci)=(-X0+idx(cj,ci));
                      vp(cj,ci)=(-Y0+idy(cj,ci));
                      xp(cj,ci)=(ii+(M/2)-1);
                      yp(cj,ci)=(jj+(N/2)-1);
                      SnR(cj,ci)=snr;
                      Pkh(cj,ci)=R(max_y1,max_x1);
                  else
                      up(cj,ci)=NaN; vp(cj,ci)=NaN; SnR(cj,ci)=NaN; Pkh(cj,ci)=0;
                      xp(cj,ci)=(ii+(M/2)-1);
                      yp(cj,ci)=(jj+(N/2)-1);
                  end
              else
                  up(cj,ci)=NaN; vp(cj,ci)=NaN; SnR(cj,ci)=NaN; Pkh(cj,ci)=0;
                  xp(cj,ci)=(ii+(M/2)-1);
                  yp(cj,ci)=(jj+(N/2)-1);
              end
              ci=ci+1;
          else
              xp(cj,ci)=(M/2)+ii-1;
              yp(cj,ci)=(N/2)+jj-1;
              up(cj,ci)=NaN; vp(cj,ci)=NaN;
              SnR(cj,ci)=NaN; Pkh(cj,ci)=NaN;ci=ci+1;
          end
          brc(cj,ci)=breakoutcounter;
      end
      breakoutcounter=1;

      % disp([num2str((cj-1)*(ci)+ci-1) ' vectors in ' num2str(toc) ' seconds'])
      fprintf('\r No. of vectors: %d', ((cj-1)*(ci)+ci-1) -sum(isnan(up(:))))
      fprintf(', Seconds taken: %f', toc);
      cj=cj+1;
  end
  fprintf('\n')
  return;
end

function [hu,hv]=globfilt(x,y,u,v,thresh)

% function [u,v]=globfilt(x,y,u,v,actions)
%
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

  fprintf(' Global filter running - ')
  norm = (sqrt(u(:).^2+v(:).^2));
  if max(norm)>0
    scale=2/max(norm);
  else
    scale=0.1;
  end

  [xo, sx] = mymean(u(:));
  [yo, sy] = mymean(v(:));

  dist = (u(:) - xo).^2 + (v(:) - yo).^2;
  valids = (dist <= (thresh*sx)^2 + (thresh*sy)^2);

  u(~valids)=NaN; v(~valids)=NaN;

  fprintf([' ..... ',num2str(sum(~valids(:))-sum(isnan(u(:)))),...
        ' vectors changed\n'])
  hu=u; hv=v;
  return;
end

function [hu,hv]=localfilt(x,y,u,v,threshold,m,maske)

% [NewU,NewV]=localfilt(x,y,u,v,threshold,method,kernelsize,mask)
%  
% This function is a filter that will remove vectors that deviate from
% the median or the mean of their surrounding neighbors by the factor
% THRESHOLD times the standard deviation of the neighbors. 
%
% METHOD (optional) should be either 'median' or 'mean'. Default is
% 'median'.
%
% KERNELSIZE is optional and is specified as number, typically 3 or 5
% which defines the number of vectors contributing to the median or mean
% value of each vector. 
%
% MASK should be applied to save calculation time. LOCALFILT is relatively 
% slow on large matrices and exlcuding "non-interesting" regions with MASK
% can increase speed in some cases.
% 
% 
% See also: matpiv, snrfilt, globfilt, peakfilt, mask


% 1999 - 2002 , jks@math.uio.no
% For use with MatPIV 1.6
%
% Copyright J.K.Sveen (jks@math.uio.no)
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
% Distributed under the Gnu General Public License
%
% Time: 10:41, Jan 17 2002

  method='mnanmedian'; stat='median'; ff=1;

  nu=NaN(size(u)+2*floor(m/2));
  nv=NaN(size(u)+2*floor(m/2));
  nu(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=u;
  nv(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=v;

  INx=false(size(nu));
  INx(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=maske;

  prev=isnan(nu); previndx=find(prev==1); 
  U2=nu+i*nv; teller=1; [ma,na]=size(U2); histo=zeros(size(nu));
  histostd=zeros(size(nu));hista=zeros(size(nu));histastd=zeros(size(nu));
  fprintf([' Local ',stat,' filter running: '])

  for ii=m-1:1:na-m+2  
      for jj=m-1:1:ma-m+2
          if INx(jj,ii)

              tmp=U2(jj-floor(m/2):jj+floor(m/2),ii-floor(m/2):ii+floor(m/2)); 
              tmp(ceil(m/2),ceil(m/2))=NaN;

              usum=nanmedian(tmp(:));

              histostd(jj,ii)=nanstd(tmp(:));
          else
              usum=nan; tmp=NaN; histostd(jj,ii)=nan;
          end
  %         u1=real(usum).^2 - real(U2(jj,ii)).^2;
  %         v1=imag(usum).^2 - imag(U2(jj,ii)).^2;
  %         
  %         histo(jj,ii)=u1+i*v1;
          histo(jj,ii)=usum;
          %histostd(jj,ii)=mnanstd(real(tmp(:))) + i*mnanstd(imag(tmp(:)));

          %th1=angle(usum); th2=angle(U2(jj,ii));
          %if th1<0, th1=2*pi+th1; end
          %if th2<0, th2=2*pi+th2; end
          %hista(jj,ii)=(th1-th2);
          %if hista(jj,ii)<0, hista(jj,ii)=2*pi+hista(jj,ii); end 
          %histastd(jj,ii)=mnanstd(abs(angle(tmp(:))));
      end
      fprintf('.')

  end

  %%%%%%%% Locate gridpoints with a higher value than the threshold 

  %[cy,cx]=find((real(histo)>threshold*real(histostd) | ...
  %    imag(histo)>threshold*imag(histostd)));

  %[cy,cx]=find( ( real(U2)>real(histo)+threshold*real(histostd) |...
  %    imag(U2)>imag(histo)+threshold*imag(histostd) |...
  %    real(U2)<real(histo)-threshold*real(histostd) |...
  %    imag(U2)<imag(histo)-threshold*imag(histostd) ) );

  bads=( real(U2)>real(histo)+threshold*real(histostd) |...
      imag(U2)>imag(histo)+threshold*imag(histostd) |...
      real(U2)<real(histo)-threshold*real(histostd) |...
      imag(U2)<imag(histo)-threshold*imag(histostd) );

  nu(bads)=NaN;
  nv(bads)=NaN;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %for jj=1:length(cy)
      %uv2(jj)=u(cy(jj),cx(jj)); vv2(jj)=v(cy(jj),cx(jj));
      %xv2(jj)=x(cy(jj),cx(jj)); yv2(jj)=y(cy(jj),cx(jj));
      % Now we asign NotANumber (NaN) to all the points in the matrix that
      % exceeds our threshold.
  %    nu(cy(jj),cx(jj))=NaN;  nv(cy(jj),cx(jj))=NaN;
  %end

  %rest=length(cy);
  rest=sum(bads(:));

  %rest2=sum(isnan(u(:)))-sum(prev(:));
  fprintf([num2str(rest),' vectors changed'])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we check for NaN's and interpolate where they exist
  %if any(strcmp(varargin,'interp'))
  %    if any(isnan(u(:)))
  %        [nu,nv]=naninterp(nu,nv);
  %    end
  %end

  hu=nu(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  hv=nv(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  fprintf('.\n')

  return;
end

function [u,v]=naninterp2(u,v,maske,xx,yy)
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

  % determine Calling m-file:
  %[stru,II]=dbstack;
  %if length(stru)>2
  %    test=stru(3).name;
  %    I2=findstr(test,'multipass');
  %    if isempty(I2), I2=0; end
  %else
  %    I2=0;
  %end 
  %if nargin==2
  %    [py,px]=find(isnan(u)==1);  
  %else
  %    py2=[];px2=[]; ipol2=zeros(size(xx));
  %    for i=1:size(mask,2)
  %        if I2~=0
  %            ipol1=inpolygon(xx,yy,mask(i).idx,mask(i).idy);
  %            ipol2=ipol2+ipol1;
  %        else
  %            ipol1=inpolygon(xx,yy,mask(i).idxw,mask(i).idyw);
  %            ipol2=ipol2+ipol1;
  %        end 
  %    end
      %[py,px]=find(isnan(u)==1 & ~ipol2 );     
      [py,px]=find(isnan(u)==1 & maske);
  %end

  numm=size(py);
  [dy,dx]=size(u);
  teller=1;
  lp=1;
  tel=1;
  % Now sort the NaN's after how many neighbors they have that are
  % physical values. Then we first interpolate those that have 8
  % neighbors, followed by 7, 6, 5, 4, 3, 2 and 1
  % use SORTROWS to sort the numbers
  fprintf(' Interpolating outliers: ')
  %pcolor(u), hold on
  while ~isempty(py)
      % check number of neighbors
      for i=1:length(py)
          %correction if vector is on edge of matrix
          corx1=0; corx2=0; cory1=0; cory2=0;
          if py(i)==1, cory1=1; cory2=0;
          elseif py(i)==dy, cory1=0; cory2=-1; end
          if px(i)==1, corx1=1; corx2=0;
          elseif px(i)==dx,  corx1=0; corx2=-1; end
          ma = u( py(i)-1+cory1:py(i)+1+cory2, px(i)-1+corx1:px(i)+1+corx2 );
          nei(i,1)=sum(~isnan(ma(:)));
          nei(i,2)=px(i);
          nei(i,3)=py(i);
      end
      % now sort the rows of NEI to interpolate the vectors with the
      % fewest spurious neighbors.
      nei=flipud(sortrows(nei,1));
      % reconstruct the sorted outlier-vectors.
      % and only interpolate the first 50% of vectors
      ind=find(nei(:,1)>=8);
      while isempty(ind)
          ind=find(nei(:,1)>=8-tel);
          tel=tel+1;
      end
      tel=1;
      py=nei(ind,3);
      px=nei(ind,2);
      for j=1:size(py,1)
          corx1=0; corx2=0; cory1=0; cory2=0;
          if py(j)==1
              cory1=1; cory2=0;
          elseif py(j)==dy
              cory1=0; cory2=-1;
          end
          if px(j)==1
              corx1=1; corx2=0;
          elseif px(j)==dx
              corx1=0; corx2=-1;
          end
          tmpu=u(py(j)-1+cory1:py(j)+1+cory2, px(j)-1+corx1:px(j)+1+corx2);
          tmpv=v(py(j)-1+cory1:py(j)+1+cory2, px(j)-1+corx1:px(j)+1+corx2);
          u(py(j),px(j))=nanmean(tmpu(:));
          v(py(j),px(j))=nanmean(tmpv(:));
          if lp>numm(1), u(py(j),px(j))=0;v(py(j),px(j))=0;end
          teller=teller+1;
      end 
      tt=length(py);

      %if nargin==2
      %    [py,px]=find(isnan(u)==1);  
      %else
          %in2=zeros(size(xx));
          %for i=1:length(mask)
          %    in=inpolygon(xx,yy,mask(i).idxw,mask(i).idyw);
          %    in2=in2+double(in);
          %end
          %[py,px]=find(isnan(u)==1 & ~ipol2 );
          [py,px]=find(isnan(u)==1 & maske);
      %end

      lp=lp+1;
      fprintf('.')
  end
  if numm(1)~=0

      fprintf([num2str(numm(1)),' Nan''s interpolated.\n'])
  else
      fprintf('Nothing to interpolate \n')
  end 

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
