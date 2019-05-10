%% code copy-pasted from  https://doi.org/10.1063/1.5016546

%%%%%%%%%%%%%%1.5 MDSR1 MATLAB Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program eliminates the strip artifacts for stripes at eliminations
%angles -90<th_elim<=90 degrees (Works only for one angle within the range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
%mex -setup;

%Make the Below path as the Current Folder
%*************************************************************
imgdir=('C:\Users\Seyyed Muhammad\Documents\MATLAB\Images4');
%*************************************************************

%maybe add test for whether imgdir exists

cd(imgdir);

%Change the format of the image if you want to e.g. "jpg"
%*************************************************************
Files = dir('*.jpg');
%*************************************************************

%Find the total number of JPEG files in the Current Folder
NumFiles= size(Files,1);

i=1;
im1 = imread(Files(i).name);
im=rgb2gray(im1);

Lx=length(im(1,:));
Ly=length(im(:,1));

% Parameters:
%*************************************************************
Ni=5;  %Number of layers(1 to 5). 5 usually gives the best results.
Nd=3;  %Power of directional decomposition, should be either "2" which creates 4 directional wedges or "3" which creates 8 directional wedges
%*************************************************************

nlevels = Nd*ones(1,Ni);        % Decomposition level e.g. [3 3 3 3] would generate 4 layers with 2^3 subimages for each
pfilter = 'maxflat'  ;          % Pyramidal filter
dfilter = 'dmaxflat7' ;         % Directional filter

% Parameters:
%*************************************************************
th_elim=(pi/180)*0;   %The elimination angle within -90<th_elim<=90 degree
S=10;  %This is "Sigma", the control over supression degree
Sa=2;  %This is "Sigma_a", the control over supression weight
%*************************************************************

th_0=abs(th_elim);
Nl=2^Nd; %Number of subiimages in each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Flipping the image for negative th_0
Ori_im=im;
if th_elim<0
im=fliplr(im);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nonsubsampled Contourlet decomposition
coeffs = nsctdec( double(im), nlevels, dfilter, pfilter );

%Adjusting the order
midpt=0.5*(Nl);
OrderVec=linspace(1,Nl,Nl);
% OrderVec = circshift(OrderVec,2^(Nd-2));
Shift=2^(Nd-2);
OrderVec = circshift(OrderVec,[1 Shift]);
OrderVec=[OrderVec(1:midpt) fliplr(OrderVec(midpt+1:Nl))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeffs_md{1,1}=coeffs{1,1};  %This is the 1st low bandpass image

for i=1:1:Ni
for j=1:1:Nl

R=coeffs{1,i+1}{1,j};
RF=fft2(R);

th(i,j)=0+(OrderVec(j)-1)*((pi)/Nl);  %This allocated the right angle to image decomposed image
delta_theta=abs(th_0-th(i,j));
SP=S*exp(-0.5*((delta_theta)^2)/(Sa^2));

for u=1:1:Ly
    for v=1:1:Lx
        Uhat(u,v)=u*cos(pi/2+th_0)+v*sin(pi/2+th_0);
        w(u,v)=1-exp(-0.5*(Uhat(u,v)^2)/(SP^2));
        RFP(u,v)=RF(u,v)*w(u,v);
    end
end

RFP=ifft2(RFP, 'symmetric');
coeffs{1,i+1}{1,j}=RFP;

end
end

%Display the coefficients and subimages
% disp('Displaying the contourlet coefficients...') ;
% shownsct( coeffs ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsubsampled Contourlet transform (NSCT) reconstruction.
% This is the inverse of nsctdec, i.e.
% imrec = nsctrec(coeffs, dfilter, pfilter);
% would reconstruct imrec = im
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reconstruct image
imrec = nsctrec( coeffs, dfilter, pfilter ) ;

%Flipping the stripe-eliminated image back if th_0 was negative
if th_elim<0
imrec=fliplr(imrec);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Displaying the reconstructed image...') ;
disp('It should be a perfect reconstruction' ) ;
disp(' ') ;

%Show the destriped/reconstructed image and the original image next to each other in
%MATLAB, you can comment it, if you don't want to see that
%*************************************************************
figure;
subplot(1,2,1), imagesc( Ori_im, [0, 255] );
title('Original image' ) ;
colormap(gray);
axis image off;
subplot(1,2,2), imagesc( imrec, [0, 255] );
title('Reconstructed image' ) ;
colormap(gray);
axis image off;
%*************************************************************

mse = sum( sum( (imrec - double(im)).^2 ) );
mse = mse / numel(im);

fprintf('The mean square error is: %f\n', mse ) ;
disp(' ');

%Saving the reconstructed image
% Normalize final
final =mat2gray(imrec,[0 255]);

% Now do the write of the destriped image, you can commnet it, if you don't
% want the destriped/reconstructed image to be saved (saves within the same folder as the
% original image)
%*************************************************************
imwrite(final,'reconstructed.png');
%*************************************************************
