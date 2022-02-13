function [ M_mean_pos ] = mean_matrix_v2( M_bound, X, Y, Z, nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


M_mean_pos = zeros(13,13,3,nbrOfFrames);

[ S1X,S2X,S3X,S4X ] = summit( X, nbrOfFrames );
[ S1Y,S2Y,S3Y,S4Y ] = summit( Y, nbrOfFrames );
[ S1Z,S2Z,S3Z,S4Z ] = summit( Z, nbrOfFrames );
[ aX,bX,cX,dX ] = plane( S1X,S2X,S3X,S4X,nbrOfFrames );
[ aY,bY,cY,dY ] = plane( S1Y,S2Y,S3Y,S4Y,nbrOfFrames );
[ aZ,bZ,cZ,dZ ] = plane( S1Z,S2Z,S3Z,S4Z,nbrOfFrames );


for k=1:nbrOfFrames
    for i=1:13
        for j=1:13
            
            left_limit  = min( M_bound(i  ,j  ,1,k) , M_bound(i+1,j  ,1,k) );
            right_limit = max( M_bound(i  ,j+1,1,k) , M_bound(i+1,j+1,1,k) );
            up_limit    = max( M_bound(i  ,j  ,2,k) , M_bound(i  ,j+1,2,k) );
            down_limit  = min( M_bound(i+1,j  ,2,k) , M_bound(i+1,j+1,2,k) );
            
            a1 = (M_bound(i+1,j  ,1,k) - M_bound(i  ,j  ,1,k)) / (M_bound(i+1,j  ,2,k) - M_bound(i  ,j  ,2,k));
            a2 = (M_bound(i+1,j+1,2,k) - M_bound(i+1,j  ,2,k)) / (M_bound(i+1,j+1,1,k) - M_bound(i+1,j  ,1,k));
            a3 = (M_bound(i  ,j+1,1,k) - M_bound(i+1,j+1,1,k)) / (M_bound(i  ,j+1,2,k) - M_bound(i+1,j+1,2,k));
            a4 = (M_bound(i  ,j  ,2,k) - M_bound(  i,j+1,2,k)) / (M_bound(i  ,j  ,1,k) - M_bound(i  ,j+1,1,k));

            b1 = M_bound(i  ,j  ,1,k) - a1 * M_bound(i  ,j  ,2,k);
            b2 = M_bound(i+1,j  ,2,k) - a2 * M_bound(i+1,j  ,1,k);
            b3 = M_bound(i+1,j+1,1,k) - a3 * M_bound(i+1,j+1,2,k);
            b4 = M_bound(  i,j+1,2,k) - a4 * M_bound(i  ,j+1,1,k);
            
            tempX = [];
            tempY = [];
            tempZ = [];
            
            
            for ix = ceil(left_limit) : floor(right_limit)
                for iy = ceil(up_limit) : floor(down_limit)
                   
                    if Z(iy,ix,k)~=0  && (ix-(a1*iy+b1))>=0 && (iy-(a2*ix+b2))<=0 && (ix-(a3*iy+b3))<=0 && (iy-(a4*ix+b4))>=0
                        
                        tempX = [tempX X(iy,ix,k)+1/cX(k)*(aX(k)*ix+bX(k)*iy+dX(k))];
                        tempY = [tempY Y(iy,ix,k)+1/cY(k)*(aY(k)*ix+bY(k)*iy+dY(k))];
                        %tempZ = [tempZ Z(iy,ix,k)+1/cZ(k)*(aZ(k)*ix+bZ(k)*iy+dZ(k))];
                        tempZ = [tempZ Z(iy,ix,k)];
                        
                    end
                    
                end
            end
            
            M_mean_pos(i,j,1,k) = mean(tempX)-1/cX(k)*(aX(k)*mean([M_bound(i,j,1,k) M_bound(i+1,j,1,k) M_bound(i+1,j+1,1,k) M_bound(i,j+1,1,k)])+bX(k)*mean([M_bound(i,j,2,k) M_bound(i+1,j,2,k) M_bound(i+1,j+1,2,k) M_bound(i,j+1,2,k)])+dX(k));
            M_mean_pos(i,j,2,k) = mean(tempY)-1/cY(k)*(aY(k)*mean([M_bound(i,j,1,k) M_bound(i+1,j,1,k) M_bound(i+1,j+1,1,k) M_bound(i,j+1,1,k)])+bY(k)*mean([M_bound(i,j,2,k) M_bound(i+1,j,2,k) M_bound(i+1,j+1,2,k) M_bound(i,j+1,2,k)])+dY(k));
            %M_mean_pos(i,j,3,k) = mean(tempZ)-1/cZ(k)*(aZ(k)*mean([M_bound(i,j,1,k) M_bound(i+1,j,1,k) M_bound(i+1,j+1,1,k) M_bound(i,j+1,1,k)])+bZ(k)*mean([M_bound(i,j,2,k) M_bound(i+1,j,2,k) M_bound(i+1,j+1,2,k) M_bound(i,j+1,2,k)])+dZ(k));
            M_mean_pos(i,j,3,k) = mean(tempZ);
            
        end
    end
    
    M_mean_pos(7,7,:,k) = nan;
    
    M_mean_pos(:,:,:,k) = 1/2*(fillmissing(M_mean_pos(:,:,:,k),'pchip',1) + fillmissing(M_mean_pos(:,:,:,k),'pchip',2));
    
end



end

