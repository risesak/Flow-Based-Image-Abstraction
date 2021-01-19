function output = Cg_one_iter(image,sigma_g,r_g,grad,m,n,p)
%% This function does smoothing along the gradient direction ,
% effectvely smoothing the interior regions of shapes.

C_g = zeros([m,n,p]);
for i=1:m
    for j=1:n
        total_wt = 0;
        
        %% This was a code block attempting speedup
%         if grad(i,j,1) > 0
%             k_min_1 = ceil((1-i)/grad(i,j,1));
%             k_max_1 = floor((m-i)/grad(i,j,1));
%         else if grad(i,j,1) < 0
%             k_max_1 = floor((1-i)/grad(i,j,1));
%             k_min_1 = ceil((m-i)/grad(i,j,1));
%             else  
%                 k_max_1 = floor(4*sigma_g);
%                 k_min_1 = -floor(4*sigma_g);
%             end        
%         end
%         
%          if grad(i,j,2) > 0
%             k_min_2 = ceil((1-j)/grad(i,j,2));
%             k_max_2 = floor((n-j)/grad(i,j,2));
%          else if grad(i,j,2) < 0
%             k_max_2 = floor((1-j)/grad(i,j,2));
%             k_min_2 = ceil((n-j)/grad(i,j,2));
%             else  
%                 k_max_2 = floor(4*sigma_g);
%                 k_min_2 = -floor(4*sigma_g);
%             end
%          end
%          
%          k_cut_neg = max(max(k_min_1,k_min_2),-floor(4*sigma_g));
%          k_cut_pos = min(min(k_max_1,k_max_2),floor(4*sigma_g));
            
        
        for k = -floor(4*sigma_g):floor(4*sigma_g)  %% need to check on the range of -T to T
            if (i+k*grad(i,j,1) < 1) || (i+k*grad(i,j,1) > m) || (j+k*grad(i,j,2) < 1) || (j+k*grad(i,j,2) > n)
                continue
            end
            dist_sqrd =  sum((k*grad(i,j,:)).^2);
            wt_spatial = exp(-dist_sqrd/(2*(sigma_g^2)));
            intensity_diff = image(i,j,:)- image(i+k*grad(i,j,1),j+k*grad(i,j,2),:);
            %% Colour dimensions handled by finding common weights for all 3 channnels with inputs from all 3 channels. 
            % Simply L2 norm is replaced in placed of modulus(for 1D case) as a distance measure for the multi channel case..
            wt_intensity = exp(-double(sum(intensity_diff.^2)/(2*(r_g^2))));
            C_g(i,j,:) = C_g(i,j,:) + wt_spatial*wt_intensity*image(i+k*grad(i,j,1),j+k*grad(i,j,2),:);
            total_wt = total_wt + wt_spatial*wt_intensity; %% Accumulate sum of weights to normalize weights later
        end
        C_g(i,j,:) = C_g(i,j,:)/total_wt;  %% Normalize weights to sum to 1
    end
end
output = C_g;                                   
end