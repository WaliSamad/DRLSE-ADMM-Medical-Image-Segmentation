%%%%these codes are modified my Dr.samad wali
clc;
close all;
clear all;
rng default
Img=imread('5bird42049.jpg');
size(Img)
    Img=rgb2gray(Img);

Img=double(Img);
%  figure;
%  imagesc(Img); colormap('gray');

[Ny,Nx] = size(Img);
% figure;imshow(Img)
%     H = fspecial('gaussian',[11,11],9);
% %    H = fspecial('average',5);
%       H = fspecial('motion',15,45);
%   sigma = 1.e-5;
%   Img = imfilter(Img,H,'circular','conv') + sigma*randn(Ny,Nx);
% Img= imnoise(Img, 'salt & pepper', 0.5); %%%worked
% Img= imnoise(Img, 'gaussian', 0,0.5);
%  figure;
%  imagesc(Img); colormap('gray');
% %  saveas(gcf, 'camel','png')


  Img = ( Img-min(Img(:)) )/ ( max(Img(:))-min(Img(:)) );



%% parameter setting
sigma=0.5;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma); % Caussian kernel
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
% region term
c1 = 0.26;
c2 = 0.91;
wrin = (Img - c2).^2;
wrout = (Img - c1).^2;
gr = wrout - wrin; % region term

size(g)
% initialize LSF as binary step function
c0=2;
initialLSF = c0*ones(size(Img));
% generate the initial region R0 as two rectangles
%%%%ground
%        initialLSF(18:28,19:24)=-c0; %%%for Gaussian Blur
  initialLSF(55:66,90:99)=-c0; %%%for Motion Blur

phi=initialLSF;

%  figure;
%  imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r'); 
%  % start level set evolution
%    axis off;
%    saveas(gcf, 'camelAdmmInitiallevelcontour','png')
alpha=0.3;
beta=4;
mu=0.1;
%%a0p3b4m0p4ru0p5rq2rv0p1rp0p01e1p1
r_u=0.5;
r_q=2;
r_varphi=0.1;
r_p=0.01;
epsilon=1.1;
phi_0=phi;
  % time step


% iter=10;%%%ok
 iter=150;

 tic
[phi,out,Dice_history] =  DRLSE_ADMM(Img,phi_0, g,gr,mu, alpha,beta, epsilon, r_u, r_q, r_varphi, r_p,iter);  
 toc
%   Dice_history
 
finalLSF=phi;
% figure;
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% hold on;  contour(phi, [0,0], 'r');
 figure;
 imagesc(Img); colormap('gray'); hold on; contour(phi,[0 0],'r'); 
 % title('Final zero level contour');
    axis off;
    
 figure;
 image(Img); colormap('white'); hold on; contour(phi,[0 0],'r'); 
  title('Final zero level contour');
axis off;
    
    





