function [tangent,grad] = ETF_to_tangent(ETF,m,n)
% calculate tangents and grad at each pixel
% tangents and grad represented using ---------------------
%                                    |(-1,-1)|(0,-1)|(1,-1)|
%                                    |(-1,0) |(0,0) |(1,0) |
%                                    |(-1,1) |(0,1) |(1,1) |
%                                     ---------------------
    tangent_pg1 = zeros([m,n]);
    tangent_pg1(-1 <= ETF(:,:,1) & ETF(:,:,1) <= -sin(pi/8)) = -1;
    tangent_pg1(-sin(pi/8) < ETF(:,:,1) & ETF(:,:,1) <= sin(pi/8)) = 0;
    tangent_pg1(sin(pi/8) < ETF(:,:,1) & ETF(:,:,1) <= 1) = 1;
    
    tangent_pg2 = zeros([m,n]);
    tangent_pg2(-1 <= ETF(:,:,2) & ETF(:,:,2) <= -sin(pi/8)) = -1;
    tangent_pg2(-sin(pi/8) < ETF(:,:,2) & ETF(:,:,2)  <= sin(pi/8)) = 0;
    tangent_pg2(sin(pi/8) < ETF(:,:,2) & ETF(:,:,2) <= 1) = 1;
    
    tangent = zeros(size(ETF));
    grad = zeros(size(ETF));
    tangent(:,:,1) = tangent_pg2;   tangent(:,:,2) = tangent_pg1;
    grad(:,:,1) = - tangent_pg1;     grad(:,:,2) = tangent_pg2;
end