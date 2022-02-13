function [M_flat] = flattener(M,a,b,c,d,nbrOfFrames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


M_flat = M;

for k=1:nbrOfFrames
    for i=1:size(M,1)
        for j=1:size(M,2)


            if numel(size(M))==3
                if M(i,j,k)~=0
                    M_flat(i,j,k) = M(i,j,k) + 1/c(k)*(a(k)*j+b(k)*i+d(k));
                end
            end

            if numel(size(M))==4
                
                n = mean([ cross(M(13,1,:,k)-M(1,1,:,k),M(1,13,:,k)-M(1,1,:,k)) , cross(M(13,13,:,k)-M(13,1,:,k),M(1,1,:,k)-M(13,1,:,k)) , cross(M(1,13,:,k)-M(13,13,:,k),M(13,1,:,k)-M(13,13,:,k)) , cross(M(1,1,:,k)-M(1,13,:,k),M(13,13,:,k)-M(1,13,:,k)) ] , 2);
    
                a(k) = n(1);
                b(k) = n(2);
                c(k) = n(3);
                d(k) = -( a(k)*mean([M(1,1,1,k),M(13,1,1,k),M(13,13,1,k),M(1,13,1,k)]) + b(k)*mean([M(1,1,2,k),M(13,1,2,k),M(13,13,2,k),M(1,13,2,k)]) + c(k)*mean([M(1,1,3,k),M(13,1,3,k),M(13,13,3,k),M(1,13,3,k)]) );
    
                if M(i,j,3,k)~=0
                    M_flat(i,j,3,k) = M(i,j,3,k) + 1/c(k)*(a(k)*M(i,j,1,k)+b(k)*M(i,j,2,k)+d(k));
                end
                
            end



        end
    end
end

end

