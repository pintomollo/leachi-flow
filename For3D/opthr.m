function T=opthr(I, print_out)
%*********************************************************************************
%FUNCTION T=OPTHR(I)
%AUTHOR: Felix Toran Marti.
%DATE: 8/5/00
%MATLAB VERSION: 5.x.
%To contact author use:
%E-MAIL: ftoran@aimme.es
%PHONE: +34 654082088
%
%This function finds the optimal threshold corresponding to the intensity image I.
%The function is intended to be a enhancement of the images toolbox for thresholding
%purposes. It can be a quick way to automate the process of manually selecting a
%threshold after seeing the histogram of an image. Also, the function helps user
%finding a reasonable good threshold value when the selection is not evident.
%
%The following example code reads a TIFF image, finds its optimal threshold, and converts
%it to a binary image:
%
%[x,map]=tiffread('c:\myimage.tiff');
%I=ind2gray(x,map);
%threshold=opthr(I);
%B=im2bw(I,threshold);
%imshow(B)
%
%If the histogram of image I is purely bimodal, the threshold will take a value
%in the middle of the valley between the 2 modes (the logical election).
%In other difficult cases, when the modes are overlapped, the threshold will minimize
%the error of interpreting background pixels as objects pixels, and vice versa.
%
%This algorithm is a small version of a more complex statistical method,
%offering good results and normally using a reduced number of iterations.
%*********************************************************************************

if nargin<2, print_out = 1; end

%Image size
[rows,cols]=size(I);

I(isnan(I))=0;

%Initial consideration: each corner of the image has background pixels.
%This provides an initial threshold (T), calculated as the mean of the gray levels contained
%in the corners. The width and height of each corner is a tenth of the image's width
%and height, respectively.


col_c=ceil(cols/10);
rows_c=ceil(rows/10);


corners=[I(1:rows_c,1:col_c); I(1:rows_c,(end-col_c+1):end);...
    I((end-rows_c+1):end,1:col_c);I((end-rows_c+1):end,(end-col_c+1):end)];


T=mean(mean(corners));


cnt = 0;
if print_out, fprintf('      '), end

%***************************************************************
% ITERATIVE PROCESS
%***************************************************************


while 1
    
    cnt = cnt+1;
    if print_out, fprintf('\b\b\b\b\b\b%6i',cnt), end
    
    %1. The mean of gray levels corresponding to objects in the image is calculated.
    %The actual threshold (T) is used to determine the boundary between objects and
    %background.
    mean_obj=sum(sum((I>T).*I ))/length(find(I>T));
    if isnan(mean_obj), mean_obj = 0; end
    
    %2. The same is done for the background pixels.
    mean_backgnd=sum(sum( (I<=T).*I ))/length(find(I<=T));
    
    %3. A new threshold is calculated as the mean of the last results:
    new_T=(mean_obj+mean_backgnd)/2;
    
    
    %4. A new iteration starts only if the threshold has changed.
    if(new_T==T)
        break;
    else
        T=new_T;
    end
    
end

if print_out, fprintf('\n'), end
%At this stage, the optimal threshold value is contained in T.