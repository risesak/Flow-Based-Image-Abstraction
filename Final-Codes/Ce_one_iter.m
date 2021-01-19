function output = Ce_one_iter(image,sigma_e,r_e,tangent,m,n,p)
%% This function does smoothing along the edge direction
C_e = zeros([m,n,p]);
for i=1:m
    for j=1:n
        % Handling every pixel separately.
        x = i; y = j;
        dist = 0;
        wt_spatial = exp(-dist.^2/(2*sigma_e*sigma_e));
        wt_intensity = exp(-(0.^2)/(2*(r_e^2)));
        % Contribution of pixel at (i,j)
        C_e(i,j,:) = C_e(i,j,:) + wt_spatial*wt_intensity*image(x,y,:);
        total_wt = wt_spatial*wt_intensity;      
        
        % Contribution of pixels along +ve edge direction
        for k = 1:floor(4*sigma_e)
            % Handling corner cases
            if x < 1 || x > m || y <1 || y > n
                break
            end
            if x + tangent(x,y,1)<1
                continue
            elseif x + tangent(x,y,1)>m
                continue
            elseif y + tangent(x,y,2)<1
                continue
            elseif y + tangent(x,y,2)>n
                continue
            end
            % Accumulate distances along the curve 
            % as direct distance would be wrong since curve may not be a straight line
            dist = dist + sqrt(tangent(x,y,1)^2 + tangent(x,y,2)^2); 
            wt_spatial = exp(-dist.^2/(2*sigma_e*sigma_e));
            intensity_diff = image(i,j,:)- image(x + tangent(x,y,1),y + tangent(x,y,2),:);
            wt_intensity = exp(-(sum(intensity_diff.^2))/(2*(r_e^2)));
            C_e(i,j,:) = C_e(i,j,:) + wt_spatial*wt_intensity*image(x + tangent(x,y,1),y + tangent(x,y,2),:);
            total_wt = total_wt + wt_spatial*wt_intensity;
            tempx = x + tangent(x,y,1);
            y = y + tangent(x,y,2);
            x = tempx;
        end
        
        % Contribution of pixels along -ve edge direction
        x = i; y = j;
        dist = 0;
        for k = 1:floor(4*sigma_e)
            % Handling corner cases
            if x < 1 || x > m || y <1 || y > n
                break
            end
            if x - tangent(x,y,1)<1
                continue
            elseif x - tangent(x,y,1)>m
                continue
            elseif y - tangent(x,y,2)<1
                continue
            elseif y - tangent(x,y,2)>n
                continue
            end            
            dist = dist + sqrt(tangent(x,y,1)^2 + tangent(x,y,2)^2);
            wt_spatial = exp(-dist.^2/(2*sigma_e*sigma_e));
            intensity_diff = image(i,j,:)- image(x - tangent(x,y,1),y - tangent(x,y,2),:);
            wt_intensity = exp(-(sum(intensity_diff.^2))/(2*(r_e^2)));
            C_e(i,j,:) = C_e(i,j,:) + wt_spatial*wt_intensity*image(x - tangent(x,y,1),y - tangent(x,y,2),:);
            total_wt = total_wt + wt_spatial*wt_intensity;
            
            tempx = x - tangent(x,y,1);
            y = y - tangent(x,y,2);
            x = tempx;
        end
        C_e(i,j,:) = C_e(i,j,:)/total_wt;  % normalize weights to sum to 1 by accumulating total_wt
        
    end
end
output = C_e;
end