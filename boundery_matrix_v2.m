function [ M_bound, E1_bound, E2_bound, E3_bound, E4_bound ] = boundery_matrix_v2( L1,L2,L3,L4,E1,E2,E3,E4,nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


M_bound = zeros(14,14,2,nbrOfFrames);

E1_bound = zeros(2,12,nbrOfFrames);
E2_bound = zeros(2,12,nbrOfFrames);
E3_bound = zeros(2,12,nbrOfFrames);
E4_bound = zeros(2,12,nbrOfFrames);

vectX=zeros(2,12,nbrOfFrames);
vectY=zeros(2,12,nbrOfFrames);

for k=1:nbrOfFrames
    
    fracL1 = [L1(k)/24 L1(k)/12*ones(1,11) L1(k)/24];
    fracL2 = [L2(k)/24 L2(k)/12*ones(1,11) L2(k)/24];
    fracL3 = [L3(k)/24 L3(k)/12*ones(1,11) L3(k)/24];
    fracL4 = [L4(k)/24 L4(k)/12*ones(1,11) L4(k)/24];
    
    E1pos = E1(:,find(E1(1,:,k)),k);
    E2pos = E2(:,find(E2(1,:,k)),k);
    E3pos = E3(:,find(E3(1,:,k)),k);
    E4pos = E4(:,find(E4(1,:,k)),k);


    % E1_bound
    last = 0;
    sum = 0;
    
    for q=1:12
        
        p = 1;
        while sum<fracL1(q)
            sum = sum + norm([ E1pos(1,last+p+1)-E1pos(1,last+p) ; E1pos(2,last+p+1)-E1pos(2,last+p) ]);
            p=p+1;
        end
        
        d = sum-fracL1(q);
        N = norm([ E1pos(1,last+p)-E1pos(1,last+p-1) ; E1pos(2,last+p)-E1pos(2,last+p-1) ]);
        a = (E1pos(1,last+p)-E1pos(1,last+p-1)) / (E1pos(2,last+p)-E1pos(2,last+p-1));
        b = E1pos(1,last+p-1) - a * E1pos(2,last+p-1);
        i = E1pos(2,last+p-1) + (N-d)/N;
        j = a * i + b;
        E1_bound(:,q,k) = [ j ; i ];

        sum = d;
        last = last+p-1;
        
    end
    
    
    % E2_bound
    last = 0;
    sum = 0;
    
    for q=1:12
        
        p = 1;
        while sum<fracL2(q)
            sum = sum + norm([ E2pos(1,last+p+1)-E2pos(1,last+p) ; E2pos(2,last+p+1)-E2pos(2,last+p) ]);
            p=p+1;
        end
        
        d = sum-fracL2(q);
        N = norm([ E2pos(1,last+p)-E2pos(1,last+p-1) ; E2pos(2,last+p)-E2pos(2,last+p-1) ]);
        a = (E2pos(2,last+p)-E2pos(2,last+p-1)) / (E2pos(1,last+p)-E2pos(1,last+p-1));
        b = E2pos(2,last+p-1) - a * E2pos(1,last+p-1);
        j = E2pos(1,last+p-1) + (N-d)/N;
        i = a * j + b;
        E2_bound(:,q,k) = [ j ; i ];

        sum = d;
        last = last+p-1;
        
    end
    
    
    % E3_bound
    last = 0;
    sum = 0;
    
    for q=1:12
        
        p = 1;
        while sum<fracL3(q)
            sum = sum + norm([ E3pos(1,last+p+1)-E3pos(1,last+p) ; E3pos(2,last+p+1)-E3pos(2,last+p) ]);
            p=p+1;
        end
        
        d = sum-fracL3(q);
        N = norm([ E3pos(1,last+p)-E3pos(1,last+p-1) ; E3pos(2,last+p)-E3pos(2,last+p-1) ]);
        a = (E3pos(1,last+p)-E3pos(1,last+p-1)) / (E3pos(2,last+p)-E3pos(2,last+p-1));
        b = E3pos(1,last+p-1) - a * E3pos(2,last+p-1);
        i = E3pos(2,last+p-1) + (N-d)/N;
        j = a * i + b;
        E3_bound(:,q,k) = [ j ; i ];

        sum = d;
        last = last+p-1;
        
    end
    
    
    % E4_bound
    last = 0;
    sum = 0;
    
    for q=1:12
        
        p = 1;
        while sum<fracL4(q)
            sum = sum + norm([ E4pos(1,last+p+1)-E4pos(1,last+p) ; E4pos(2,last+p+1)-E4pos(2,last+p) ]);
            p=p+1;
        end
        
        d = sum-fracL4(q);
        N = norm([ E4pos(1,last+p)-E4pos(1,last+p-1) ; E4pos(2,last+p)-E4pos(2,last+p-1) ]);
        a = (E4pos(2,last+p)-E4pos(2,last+p-1)) / (E4pos(1,last+p)-E4pos(1,last+p-1));
        b = E4pos(2,last+p-1) - a * E4pos(1,last+p-1);
        j = E4pos(1,last+p-1) + (N-d)/N;
        i = a * j + b;
        E4_bound(:,q,k) = [ j ; i ];

        sum = d;
        last = last+p-1;
        
    end
    
    
    
    
    %VectXY
    for p=1:12
        vectX(:,p,k) = [(E3_bound(1,p,k)-E1_bound(1,p,k))/12 ; (E3_bound(2,p,k)-E1_bound(2,p,k))/12];
        vectY(:,p,k) = [(E2_bound(1,p,k)-E4_bound(1,p,k))/12 ; (E2_bound(2,p,k)-E4_bound(2,p,k))/12];
    end
    
    
    
    
    %M_bound
    
    M_bound( 1, 1,:,k) = [E1pos(1,  1)-vectX(1, 1,k)/2-vectY(1, 1,k)/2 ; E1pos(2,  1)-vectX(2, 1,k)/2-vectY(2, 1,k)/2];
    M_bound(14, 1,:,k) = [E2pos(1,  1)-vectX(1,12,k)/2+vectY(1, 1,k)/2 ; E2pos(2,  1)-vectX(2,12,k)/2+vectY(2, 1,k)/2];
    M_bound(14,14,:,k) = [E3pos(1,end)+vectX(1,12,k)/2+vectY(1,12,k)/2 ; E3pos(2,end)+vectX(2,12,k)/2+vectY(2,12,k)/2];
    M_bound( 1,14,:,k) = [E4pos(1,end)+vectX(1, 1,k)/2-vectY(1,12,k)/2 ; E4pos(2,end)+vectX(2, 1,k)/2-vectY(2,12,k)/2];
    
    for p=2:13
        M_bound( p, 1,:,k) = E1_bound(:,p-1,k) - vectX(:,p-1,k)/2;
        M_bound(14, p,:,k) = E2_bound(:,p-1,k) + vectY(:,p-1,k)/2;
        M_bound( p,14,:,k) = E3_bound(:,p-1,k) + vectX(:,p-1,k)/2;
        M_bound( 1, p,:,k) = E4_bound(:,p-1,k) - vectY(:,p-1,k)/2;
    end
     
%     for i=2:13
%         for j=2:13
%             M_bound(i,j,1,k) = 1/2 * (M_bound(i,1,1,k) + (j-1)*vectX(1,i-1,k) + M_bound(1,j,1,k) + (i-1)*vectY(1,j-1,k));
%             M_bound(i,j,2,k) = 1/2 * (M_bound(i,1,2,k) + (j-1)*vectX(2,i-1,k) + M_bound(1,j,2,k) + (i-1)*vectY(2,j-1,k));
%         end
%     end
    
    for i=2:13
        for j=2:13
            M_bound(i,j,1,k) = M_bound(i,1,1,k) + (j-1)*vectX(1,i-1,k);
            M_bound(i,j,2,k) = M_bound(1,j,2,k) + (i-1)*vectY(2,j-1,k);
        end
    end
    
    
end

end

