function output_image = myFBLfilter(image,ETF,sigma_g, r_g, sigma_e,r_e,num_iter)
% sigma_g : gradient direction,spatial guassian variance
% r_g : gradient direction,intensity based guassian variance
% sigma_e : edge direction,spatial guassian variance
% r_e : edge direction,intensity based guassian variance
% num_iter : num of iterations where one iteration opearates 
%both along edge and gradient once each. 
[m, n, p] = size(image);
[tangent,grad] = ETF_to_tangent(ETF,m,n);
outCe = image;
%% Alternate Smoothing along edge flow and along gradient direction
for i = 1:num_iter
    % calculate C_g at each pixel
    outCg = Cg_one_iter(outCe,sigma_g,r_g,grad,m,n,p);  % Along graident
    disp("Cg done")
    % calculate C_e at each pixel
    outCe = Ce_one_iter(outCg,sigma_e,r_e,tangent,m,n,p); % Along edge
    disp("Ce done")
end
output_image = outCe;
end





