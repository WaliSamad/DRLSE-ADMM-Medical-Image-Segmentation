%%%%%these codes are modified my Dr.samad wali

clc;
close all;
clear all;
 rng default
Img=imread('bushes.jpg');
%    Img=rgb2gray(Img);
Img=double(Img);
[Ny,Nx] = size(Img);
% figure;imshow(Img)
     Img = ( Img-min(Img(:)) )/ ( max(Img(:))-min(Img(:)) );

%% parameter setting
sigma=0.02;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma); % Caussian kernel
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
size(g)
c1 = 0.26;
c2 = 0.91;
wrin = (Img - c2).^2;
wrout = (Img - c1).^2;
gr = wrout - wrin; % region term
% % initialize LSF as binary step function
c0=2;
initialLSF = c0*ones(size(Img));
% generate the initial region R0 as two rectangles
%%%%%%twocells
initialLSF(11:26, 29:40)=-c0;

% initialLSF(21:28, 54:59)=-c0;
phi=initialLSF;

figure;
 imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r');
 title('Initial zero level contour');
 axis off;
%      saveas(gcf, 'Bushesinitiallevelsetcontour','png')
alpha=0.9;
beta=-1.5;
mu=1;

r_u=0.5;
r_q=3;
r_varphi=0.1;
r_p=0.5;
epsilon=1;
phi_0=phi;
  % time step


% iter_inner=10;%%%ok
%iter_inner=100;
iter=100;
%  wr=g;
% tic
% [phi,out,Dice_history,Jaccard_history] =  DRLSE_ADMM(Mask,Img,phi_0, g,gr,mu, alpha,beta, epsilon, r_u, r_q, r_varphi, r_p,iter_inner);  
%  toc
%   tic
[phi,out,Dice_history] =  DRLSE_ADMM(Img,phi_0, g',gr,mu, alpha,beta, epsilon, r_u, r_q, r_varphi, r_p,iter);  
 toc
  
%   Jaccard_history(end)
finalLSF=phi;
% figure;
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% hold on;  contour(phi, [0,0], 'r');
% title('Final zero level contour');
 figure;
imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r'); 
 title('Final zero level contour');
      axis off;
