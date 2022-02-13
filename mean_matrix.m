function [ M_mean_pos ] = mean_matrix( M_bound, X, Y, Z, nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


M_mean_pos = zeros(13,13,3,nbrOfFrames);

for k=1:nbrOfFrames
    for i=1:13
        for j=1:13
            
            left_limit  = max( M_bound(i  ,j  ,1,k) , M_bound(i+1,j  ,1,k) );
            right_limit = min( M_bound(i  ,j+1,1,k) , M_bound(i+1,j+1,1,k) );
            up_limit    = min( M_bound(i  ,j  ,2,k) , M_bound(i  ,j+1,2,k) );
            down_limit  = max( M_bound(i+1,j  ,2,k) , M_bound(i+1,j+1,2,k) );
            
            tempX = [];
            tempY = [];
            tempZ = [];
            
            for ix = ceil(left_limit) : floor(right_limit)
                for iy = ceil(up_limit) : floor(down_limit)
                    
                    if Z(iy,ix,k)~=0
                        
                        tempX = [tempX X(iy,ix,k)];
                        tempY = [tempY Y(iy,ix,k)];
                        tempZ = [tempZ Z(iy,ix,k)];
                        
                    end
                end
            end
            
            M_mean_pos(i,j,1,k) = mean(tempX);
            M_mean_pos(i,j,2,k) = mean(tempY);
            M_mean_pos(i,j,3,k) = mean(tempZ);
            
        end
    end
    
    M_mean_pos(7,7,:,k) = nan;
    
end



end

