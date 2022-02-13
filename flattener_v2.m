function [M_flat] = flattener_v2(M,a,b,c,d,nbrOfFrames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


M_flat = M;
M_center = M;
        

for k=1:nbrOfFrames
    
    if numel(size(M))==3
        for i=1:size(M,1)
            for j=1:size(M,2)

                if M(i,j,k)~=0
                    M_flat(i,j,k) = M(i,j,k) + 1/c(k)*(a(k)*j+b(k)*i+d(k));
                end
            end
        end
    end
            
            
    if numel(size(M))==4
        
        %Centering
        center = mean([squeeze(M(1,1,:,k))';squeeze(M(end,1,:,k))';squeeze(M(end,end,:,k))';squeeze(M(1,end,:,k))']);
        
        for i=1:size(M,1)
            for j=1:size(M,2)
                M_center(i,j,:,k) = squeeze(M(i,j,:,k))-center';
            end
        end
        
        
        n = squeeze(mean([ cross(M_center(end,1,:,k)-M_center(1,1,:,k),M_center(1,end,:,k)-M_center(1,1,:,k)) , cross(M_center(end,end,:,k)-M_center(end,1,:,k),M_center(1,1,:,k)-M_center(end,1,:,k)) , cross(M_center(1,end,:,k)-M_center(end,end,:,k),M_center(end,1,:,k)-M_center(end,end,:,k)) , cross(M_center(1,1,:,k)-M_center(1,end,:,k),M_center(end,end,:,k)-M_center(1,end,:,k)) ] , 2));
        n = n/norm(n);
        if n(3)<0
            n = -n;
        end
        
        %%
        
%                 a(k) = n(1);
%                 b(k) = n(2);
%                 c(k) = n(3);
%                 d(k) = -( a(k)*mean([M_center(1,1,1,k),M_center(end,1,1,k),M_center(end,end,1,k),M_center(1,end,1,k)]) + b(k)*mean([M_center(1,1,2,k),M_center(end,1,2,k),M_center(end,end,2,k),M_center(1,end,2,k)]) + c(k)*mean([M_center(1,1,3,k),M_center(end,1,3,k),M_center(end,end,3,k),M_center(1,end,3,k)]) );
%     
%                 for i=1:size(M,1)
%                     for j=1:size(M,2)
% %                         if M_center(i,j,3,k)~=0
%                             M_center(i,j,3,k) = M_center(i,j,3,k) + 1/c(k)*(a(k)*M_center(i,j,1,k)+b(k)*M_center(i,j,2,k)+d(k));
% %                         end
%                     end
%                 end

        

        %%
        u = cross(n,[0 0 1]);
        u = u/norm(u);
        theta = acos(n(3));

        P = u*u';
        I = eye(3);
        Q = [ 0 -u(3) u(2) ; u(3) 0 -u(1) ; -u(2) u(1) 0 ];

        R = P + cos(theta)*(I-P) + sin(theta)*Q;
        %R = eye(3);
        
        for i=1:size(M,1)
            for j=1:size(M,2)
                M_center(i,j,:,k) = real(R*squeeze(M_center(i,j,:,k)));
            end
        end
        
    
        %%
        M_flat=M_center;
        
    
    
    end
            
end

end

