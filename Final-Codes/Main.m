%% MyMainScript
tic;
clear all;
close all;

%% Reading Images
image_1 = double(imread('../data/forest.jpg'));
% image_1 = image_1(1:2:end,1:2:end,:);
image_1 = image_1(250:450,20:280,:);
%% Tuning the images
%% Image1
image = rgb2gray(uint8(image_1));
figure
imshow(uint8(image_1),[]);
title('Original Image');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_size = 5;
[ETF] = ETF(image, window_size,2);
disp("ETF done")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_g = 0.3;
r_g = 10;
sigma_e = 2;
r_e = 50;
num_iter = 5;
smoothened_1 = myFBLfilter(image_1,ETF,sigma_g, r_g, sigma_e,r_e,num_iter);

figure
imshow(uint8(smoothened_1),[]);
title('smoothening(FBL filter)');
colorbar
% 
% sigma_spatial = 3;
% sigma_intensity = 22 ;
% window_size = 9;
% smoothened_2 = image_1;
% for i=1:5
% smoothened_2 = BilateralFilter(smoothened_2,sigma_spatial,sigma_intensity,window_size);
% disp("iter1 done")
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% imshow(uint8(smoothened_2),[]);
% title('smoothening(Normal Bilateral Filter)');
% colorbar 

toc;

% sigma_c = 1;
% rho = 0.99;
% sigma_m =  1;
% tau = 0.5;
% num_iter = 3;
% sigma_blur = 2;
% [H_e] = FDOG(image,ETF,sigma_c,rho,sigma_m,tau,sigma_blur,num_iter);
% figure
% imshow(uint8(H_e),[]);
% title('FDOG out');
% colorbar

% figure
% imshow(uint8(smoothened_1.*repmat(H_e,1,1,3)),[]);
% colormap jet
% title('Integrated');
% % colorbar  
