function [output_image]  = BilateralFilter(image, sigma_spatial, sigma_intensity, window_size)

[x, y, p] = size(image);
w_tra = round((window_size-1)/2);
%Now we want to use 2 for loops, so in order to be able to generalize the
%process over the entire image, we shall pad the boundaries with border elements of the original 
%array according to the value of window_size so that generalized method can be used
%for the pixels at the edges as well
img_padded = padarray(image, [w_tra,w_tra], 'replicate');
%The pixels of concern in corr_img_padded are the ones having x coordinates
%ranging from w_tra+1 to w_tra+x and y coordinates ranging from w_tra+1 to w_tra+y. where x and y
%are the size dimensions of the original input image.
% So we shall use two for loops to traverse over each pixel of the image
output_image = zeros(x,y,p);
for i = w_tra+1:w_tra+x    
    for j = w_tra+1:w_tra+y    
        % Now, for every value of i and j, we need to consider a 
        % window_size X window_size large region and calculate weights for both the
        % spatial gaussian and the intensity gaussian
    
        % For calculating weights corresponding to the intensity gaussian, we
        % will look at intensity values at pixel(i,j) with respect to the
        % intensity values at all other locations in the window_size X
        % window_size region
    
        % For calculating weights corresponding to spatial gaussian, we shall
        % look at the spatial distance of the location (i,j) with respect to
        % all the other pixel locations in the window_size X
        % window_size region
        
        % Now we are given in the question that we are not allowed to use any
        % more for loops. So we need to calculate the gaussian weights for the
        % entire window at once.
    
        % So now we define two variables which for a given (i,j) will span the
        % window_size X window_size region at once
    
        win_x = i-w_tra:i+w_tra;
        win_y = j-w_tra:j+w_tra;
    
        % Handling colour with L2 norm as the distance measure. 
        % Common weights for all 3 channels calculated from inputs from all 3 channels.        
        
        weights_intensity = exp(-sum((img_padded(win_x, win_y,:) - img_padded(i, j,:)).^2,3)/(2*(sigma_intensity)^2));
        weights_intensity = weights_intensity/sqrt(2*pi*sigma_intensity^2)  ; 
          
        %% For weights based on spatial distance, form matrices X and Y of 
        % size = (window_size,window_size) which store x and y 
        %distance of each point in the window from the centre point. 
        X =(1:window_size)'*ones(1,window_size);
        Y = ones(window_size,1)*(1:window_size);
        % shifting to make our centre point (0,0)
        X = X - (w_tra+1);  
        Y = Y - (w_tra+1);
        weights_spatial = exp(-(X.^2 + Y.^2)/(2*(sigma_spatial)^2));
        weights_spatial = weights_spatial/sqrt(2*pi*sigma_spatial^2);        
        wts = weights_intensity .* weights_spatial;
        w_p = sum(wts, "all");  % Normaliation term calculated
%         disp(w_p)
        I = img_padded(win_x, win_y,:) .* repmat(weights_intensity,1,1,p) .*  repmat(weights_spatial,1,1,p);
%         disp(repmat(weights_spatial,1,1,p))
        output_image(i-w_tra, j-w_tra,:) = sum(I,[1,2])./w_p; % sum only along 2 dimesnsions i.e. window width and height dimesnsions and keep colour dimesnsion intact.  
    end    
end
end
