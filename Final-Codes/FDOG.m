function [H_e] = FDOG(image,ETF_image,sigma_c,rho,sigma_m,tau,sigma_blur,num_iter)
%FDOG Implementation of the FDOG filter
%Inputs:
%   image: an mxn grayscale image
%   ETF_image: an mxnx2 image where,
%       ETF[.][.][1]: x-coordinate(horizontal)
%       ETF[.][.][2]: y-coordinate(vertical)
%   sigma_c: same as var used in paper for DOG along gradient
%   rho: same as var used in paper for DOG along gradient
%   sigma_m: used to weigh pixels along tangent flow
%   tau: used for thresholding
%   sigma_blur: used to blur input image before every call to FDOG_one_iter
%   num_iter: number of iterations of the FDOG filter
%Outputs:
%   H_e: output image showing edges
image = double(image);
[m,n,~] = size(ETF_image);

% unit normalise ETF
ETF_mag = sqrt(ETF_image(:,:,1).^2 + ETF_image(:,:,2).^2);
ETF_mag(ETF_mag==0) = 1;
ETF_image = ETF_image./ETF_mag;

[tangent,grad] = ETF_to_tangent(ETF_image,m,n); % get tangents and gradient neighbours from ETF

filterSize_blur = ceil(3*sigma_blur+4);
guass1= fspecial("gaussian",filterSize_blur,sigma_blur);  %compute gaussian filter

image = imfilter(image,guass1);  % blurring the image
H_e = FDOG_one_iter(image,tangent,grad,sigma_c,rho,sigma_m,tau); % one iteration of FDOG
figure
imshow(uint8(H_e),[]);
title('FDOG out');
colorbar
for i=1:num_iter-1
    H_e = H_e.*image; % add edges to image
    H_e = imfilter(H_e,guass1);  % blurring the image
    H_e = FDOG_one_iter(H_e,tangent,grad,sigma_c,rho,sigma_m,tau); % one iteration of FDOG
    figure
    imshow(uint8(H_e),[]);
    title('FDOG out');
    colorbar
end

end

function [tangent,grad] = ETF_to_tangent(ETF_image,m,n)
% calculate tangents and grad at each pixel
% tangents and grad represented using ---------------------
%                                    |(-1,-1)|(0,-1)|(1,-1)|
%                                    |(-1,0) |(0,0) |(1,0) |
%                                    |(-1,1) |(0,1) |(1,1) |
%                                     ---------------------    
    tangent_pg1 = zeros([m,n]);
    tangent_pg1((-1 <= ETF_image(:,:,1)) & (ETF_image(:,:,1) <= -sin(pi/8))) = -1;
    tangent_pg1((-sin(pi/8) < ETF_image(:,:,1)) & (ETF_image(:,:,1) <= sin(pi/8))) = 0;
    tangent_pg1((sin(pi/8) < ETF_image(:,:,1)) & (ETF_image(:,:,1) <= 1)) = 1;
    
    tangent_pg2 = zeros([m,n]);
    tangent_pg2((-1 <= ETF_image(:,:,2)) & (ETF_image(:,:,2) <= -sin(pi/8))) = -1;
    tangent_pg2((-sin(pi/8) < ETF_image(:,:,2)) & (ETF_image(:,:,2) <= sin(pi/8))) = 0;
    tangent_pg2((sin(pi/8) < ETF_image(:,:,2)) & (ETF_image(:,:,2) <= 1)) = 1;
    
    tangent = zeros(size(ETF_image));
    grad = zeros(size(ETF_image));
    tangent(:,:,1) = tangent_pg1;   tangent(:,:,2) = tangent_pg2;
    grad(:,:,1) = tangent_pg2;     grad(:,:,2) = -tangent_pg1;
end