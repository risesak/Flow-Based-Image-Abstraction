function [ETF_output]  = ETF(image,window_size,num_iter)

image = double(image);
maxi = max(max(image));
image = image/maxi;

%applying gaussian blurring to reduce noise
sigma_blur = 2;
filterSize_blur = ceil(3*sigma_blur+4);
guass1= fspecial("gaussian",filterSize_blur,sigma_blur);  %compute gaussian filter
blurred_image = imfilter(image,guass1);  % blurring the image

% Gradient computation using sobel filter
[Gx,Gy] = imgradientxy(blurred_image);
G = sqrt(Gx.^2 + Gy.^2);
G_norm = G/max(max(G));

G(G==0) = 1;

% vectors perpendicular to gradient(local tangential direction)
tx = (-1 * Gy)./G;
ty = Gx./G;

% stacking these two(tx and ty) gives the tangent vector at that particular pixel

% refer ETF_one_iter to see the algorithm for updation of tx and ty
ETF_output = ETF_one_iter(tx,ty,G_norm,window_size);
for i=1:num_iter-1
    ETF_output = ETF_one_iter(ETF_output(:,:,1),ETF_output(:,:,2),G_norm,window_size);
end

figure
imshow(uint8(ETF_output(:,:,1)),[]);
title('ETF channel 1');
colorbar
figure
imshow(uint8(ETF_output(:,:,2)),[]);
title('ETF channel 2');
colorbar


end
