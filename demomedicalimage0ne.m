%%%%these codes are modified my Dr.samad wali

clc;
close all;
clear all;
rng default
Img=imread('0104.bmp');
%   Img=rgb2gray(Img);
Img=double(Img);
[Ny,Nx] = size(Img)
% figure;imshow(Img)
  Img = ( Img-min(Img(:)) )/ ( max(Img(:))-min(Img(:)) );
% region term
c1 = 0.26;
c2 = 0.91;
wrin = (Img - c2).^2;
wrout = (Img - c1).^2;
gr = wrout - wrin; % region term


sigma=0.2;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma); % Caussian kernel
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
size(g)
% initialize LSF as binary step function
c0=2;
initialLSF = c0*ones(size(Img));
% generate the initial region R0 as two rectangles

  initialLSF(20:27, 22:26)=-c0;
phi=initialLSF;
 figure;
 imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r'); 
 title('Final zero level contour');
     axis off;
%      saveas(gcf, 'MedicalimageoneAdmmInitiallevelcontour','png')
alpha=0.9;
beta=4;
mu=1;

r_u=0.1;
r_q=3;
r_varphi=0.19;
r_p=3;
epsilon=1;
phi_0=phi;
  % time step


% iter=15;%%%ok
iter=100;

 tic
[phi,out,Dice_history] =  DRLSE_ADMM(Img,phi_0, g,gr,mu, alpha,beta, epsilon, r_u, r_q, r_varphi, r_p,iter);  
 toc
finalLSF=phi;
% figure;
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% hold on;  contour(phi, [0,0], 'r');
 figure;
 imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r'); 
% title('Final zero level contour');
  axis off;
%   saveas(gcf, 'MedicalimageoneAdmmFinallevelcontour','png')

  figure;
 image(Img); colormap('white'); hold on; contour(phi,[0 0],'r'); 
 title('Final zero level contour');
   axis off;




