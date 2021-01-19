function out = DOG(image,rho,sigma_blur)
image = double(image);
%% Find Guassian 1
filterSize_blur = ceil(3*sigma_blur+4);
guass1= fspecial("gaussian",filterSize_blur,sigma_blur);  %compute gaussian filter
image1 = imfilter(image,guass1);  % blurring the image

%% Find Guassian 2
sigma_blur = 1.6*sigma_blur;  % To approximately make it LOG
filterSize_blur = ceil(3*sigma_blur+4);
guass2= fspecial("gaussian",filterSize_blur,sigma_blur);  %compute gaussian filte
image2 = imfilter(image,guass2);  % blurring the image

%% Find DOG
out = image1 - rho*image2;
thresh = 0;

%% Used threshold's similar to the paper suggested for FDOG
out(out <= thresh & 1+ tanh(out) < 0.5) = 0;
out(out > thresh) = 1;
end

