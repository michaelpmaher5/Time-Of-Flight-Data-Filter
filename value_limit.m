function [ index,limit ] = value_limit( Z,a,b,c,d,nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


index_matrix = zeros(120,160,nbrOfFrames);

limit = zeros(nbrOfFrames,1);


for k=1:nbrOfFrames
    
    [~,~,Z_temp] = find(Z(:,:,k));
    limit(k) = 0.05*mean(Z_temp);
    %limit(k) = 2*std(Z(find(Z(:,:,k))),1);
    
    for i=1:120
        for j=1:160
            
            if Z(i,j,k)~=0 && ( Z(i,j,k)<-1/c(k)*(a(k)*j+b(k)*i+d(k))-limit(k) || Z(i,j,k)>-1/c(k)*(a(k)*j+b(k)*i+d(k))+limit(k) )
                index_matrix(i,j,k) = 1;
            end
            
        end
    end
    
end

index = (index_matrix==1);
    
end

