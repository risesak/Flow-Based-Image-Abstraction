function [H_e] = FDOG_one_iter(image,tangent,grad,sigma_c,rho,sigma_m,tau)
%FDOG_ONE_ITER Does one iteration of the FDOG filter
%Inputs:
%   image: an mxn grayscale image
%   tangent: an mxnx2 image, showing neighbor along edge tangent at each pixel, where,
%       tangent[.][.][1]: x-coordinate(horizontal)
%       tangent[.][.][2]: y-coordinate(vertical)
%   grad: an mxnx2 image, showing neighbor along grad at each pixel, where,
%       grad[.][.][1]: x-coordinate(horizontal)
%       grad[.][.][2]: y-coordinate(vertical)
%   sigma_c: same as var used in paper for DOG along gradient
%   rho: same as var used in paper for DOG along gradient
%   sigma_m: used to weigh pixels along tangent flow
%   tau: used for thresholding
%Outputs:
%   H_e: output image showing edges
[m,n,~] = size(tangent);

% calculate H_g at each pixel
H_g = zeros([m,n]);
total_wt = 0; % for use at each pixel
sigma_s = 1.6*sigma_c;
for i=1:m
    for j=1:n
        for k = -floor(4*sigma_s):floor(4*sigma_s)
            % x = j+k*grad(i,j,1); y = i+k*grad(i,j,2)
            if (j+k*grad(i,j,1)) < 1
                continue
            elseif (j+k*grad(i,j,1)) > n
                continue
            elseif (i+k*grad(i,j,2)) < 1
                continue
            elseif (i+k*grad(i,j,2)) > m
                continue
            end
            % now 1<=x<=n, 1<=y<=m
            dist_sqrd = (k*grad(i,j,1))^2 + (k*grad(i,j,2))^2;
            wt = exp(-dist_sqrd/(2*(sigma_c^2)))/sigma_c - rho*exp(-dist_sqrd/(2*(sigma_s^2)))/sigma_s;
            wt = wt/sqrt(2*pi);
            H_g(i,j) = H_g(i,j) + wt*image(i+k*grad(i,j,2),j+k*grad(i,j,1));
            total_wt = total_wt+wt;
        end
        H_g(i,j) = H_g(i,j)/total_wt;
        total_wt = 0;
    end
end

% calculate H_e at each pixel and threshold
H_e = zeros([m,n]);
for i=1:m
    for j=1:n
        x = j; y = i;
        dist = 0;
        wt = exp(-(dist^2)/(2*sigma_m*sigma_m))/(sigma_m*sqrt(2*pi));
        H_e(i,j) = H_e(i,j) + wt*H_g(y,x);
        total_wt = wt;
        for k = 1:floor(4*sigma_m)
            if (x + tangent(y,x,1))<1
                continue
            elseif (x + tangent(y,x,1))>n
                continue
            elseif (y + tangent(y,x,2))<1
                continue
            elseif (y + tangent(y,x,2))>m
                continue
            end
            % now 1<= x + tangent(y,x,1) <=n and 1<= y + tangent(y,x,2) <=m
            dist = dist + sqrt(tangent(y,x,1)^2 + tangent(y,x,2)^2);
            wt = exp(-(dist^2)/(2*sigma_m*sigma_m))/(sigma_m*sqrt(2*pi));
            H_e(i,j) = H_e(i,j) + wt*H_g(y + tangent(y,x,2),x + tangent(y,x,1));
            total_wt = total_wt + wt;
            temp_x = x; temp_y=y;
            x = temp_x + tangent(temp_y,temp_x,1);
            y = temp_y + tangent(temp_y,temp_x,2);
        end
        
        x = j; y = i;
        dist = 0;
        for k = 1:floor(4*sigma_m)
            if (x - tangent(y,x,1))<1
                continue
            elseif (x - tangent(y,x,1))>n
                continue
            elseif (y - tangent(y,x,2))<1
                continue
            elseif (y - tangent(y,x,2))>m
                continue
            end
            % now 1<= x - tangent(y,x,1) <=n and 1<= y - tangent(y,x,2) <=m
            dist = dist + sqrt(tangent(y,x,1)^2 + tangent(y,x,2)^2);
            wt = exp(-(dist^2)/(2*sigma_m*sigma_m))/(sigma_m*sqrt(2*pi));
            H_e(i,j) = H_e(i,j) + wt*H_g(y - tangent(y,x,2),x - tangent(y,x,1));
            total_wt = total_wt + wt;
            temp_x = x; temp_y=y;
            x = temp_x - tangent(temp_y,temp_x,1);
            y = temp_y - tangent(temp_y,temp_x,2);
        end
        H_e(i,j) = H_e(i,j)/total_wt;
        % threshold now
        if (H_e(i,j)<0) && (1+tanh(H_e(i,j))<tau)
            H_e(i,j) = 0;
        else
            H_e(i,j) = 1;
        end
    end
end

end
