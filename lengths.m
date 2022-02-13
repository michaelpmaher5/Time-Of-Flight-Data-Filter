function [ L1,L2,L3,L4 ] = lengths( E1,E2,E3,E4,nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


L1 = zeros(1,nbrOfFrames);
L2 = zeros(1,nbrOfFrames);
L3 = zeros(1,nbrOfFrames);
L4 = zeros(1,nbrOfFrames);


for k=1:nbrOfFrames
    
    for p=find(E1(1,:,k),numel(find(E1(1,:,k)))-1)
        L1(k) = L1(k) + norm([ E1(1,p+1,k)-E1(1,p,k) ; E1(2,p+1,k)-E1(2,p,k) ]);
    end
    
    for p=find(E2(1,:,k),numel(find(E2(1,:,k)))-1)
        L2(k) = L2(k) + norm([ E2(1,p+1,k)-E2(1,p,k) ; E2(2,p+1,k)-E2(2,p,k) ]);
    end
    
    for p=find(E3(1,:,k),numel(find(E3(1,:,k)))-1)
        L3(k) = L3(k) + norm([ E3(1,p+1,k)-E3(1,p,k) ; E3(2,p+1,k)-E3(2,p,k) ]);
    end
    
    for p=find(E4(1,:,k),numel(find(E4(1,:,k)))-1)
        L4(k) = L4(k) + norm([ E4(1,p+1,k)-E4(1,p,k) ; E4(2,p+1,k)-E4(2,p,k) ]);
    end

end

