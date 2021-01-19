function [ETF_output2]  = ETF_one_iter(tx,ty,G_norm,window_size)
%...

[x,y,~] = size(tx);
w_tra = round((window_size-1)/2);

%padding so as to avoid handling the corner courses separately
G_padded = padarray(G_norm, [w_tra,w_tra], 'replicate');
tx_padded = padarray(tx, [w_tra,w_tra], 'replicate');
ty_padded = padarray(ty, [w_tra,w_tra], 'replicate');

%initalizing array for storing output
output_image_x = zeros(x,y);
output_image_y = zeros(x,y);

%looping over allthe pixels of the image
for i = w_tra+1:w_tra+x
    
    for j = w_tra+1:w_tra+y
        
        % defining a window around the center pixel
        win_x = i-w_tra:i+w_tra;
        win_y = j-w_tra:j+w_tra;
        
        % Here G_padded(x,y) denotes the magnitude of gradient vector at
        % point (x,y). So, the weight wm is  monotonically increases with respect to the magnitude
        % difference of gradient at center pixel wrt neighboruing pixel, indicating that bigger weights
        % are given to the neighboring pixels y whose gradient magnitudes are higher than that
        % of the center x. 
        % This ensures the "preservation of the dominant edge directions"
        wm = (G_padded(win_x,win_y)- G_padded(i,j) + 1)/2;
        
        
        % X_tx and Y_tx represent the x and y components of the tangent vector
        % for the center pixel
        X_tx = tx_padded(i,j);
        Y_tx = ty_padded(i,j);
        
        % X_tx and Y_tx represent the x and y components of the tangent vector
        % for the neighbouring pixels
        X_ty = tx_padded(win_x,win_y);
        Y_ty = ty_padded(win_x,win_y);
        
        % Calculating dot product between the two tangent vectors described
        % above
        prod_x = X_tx*X_ty;
        prod_y = Y_tx*Y_ty;
        dot_prod = prod_x + prod_y;
        
        % wd is then assigned the absolute value of the dot product
        % calculated above. This weight function increases as the 
        % two vectors align closely. SO this gives more weights to those edges in neighbourhood
        % that are aligned more similar to the center pixel thereby smoothing among the edges 
        % with similar orientations.
        wd = abs(dot_prod);
        
        % For tight alignment of vectors, we reverse the
        % direction of tangent vector in the neighbour by multiplying by
        % phi, where phi takes value -1 when dot product is negative and 1
        % otherwise
        phi = dot_prod;
        phi(phi<=0)=-1;
        phi(phi>0)=1;
        
        % Multiplying all the weight matrices for the window around the
        % center pixel for both x and y component
        Ix = wm.*wd.*phi.*X_ty;
        Iy = wm.*wd.*phi.*Y_ty;
        
        output_image_x(i-w_tra, j-w_tra) = (sum(Ix, 'all'));
        output_image_y(i-w_tra, j-w_tra) = (sum(Iy, 'all'));     
    end    
end
t = sqrt(output_image_x.^2 + output_image_y.^2);
% normalizing the outputs
output_image_x = output_image_x./t;
output_image_y = output_image_y./t;
ETF_output2 = zeros(x,y,2);
% final output for this iteration
ETF_output2(:,:,1) = output_image_x;
ETF_output2(:,:,2) = output_image_y;
end
