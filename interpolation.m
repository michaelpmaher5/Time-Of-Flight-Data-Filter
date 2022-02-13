function [ Mq ] = interpolation( M,sizeX,sizeY,nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[X,Y] = ndgrid(linspace(1,size(M,1),sizeX),linspace(1,size(M,2),sizeY));

Mq = zeros(size(X,1),size(Y,2),3,nbrOfFrames);


for k=1:nbrOfFrames

    F = griddedInterpolant(M(:,:,3,k),'makima');
    Zq = F(X,Y);

    F.Values = M(:,:,1,k);
    Xq = F(X,Y);
    F.Values = M(:,:,2,k);
    Yq = F(X,Y);

    
    Mq(:,:,1,k) = Xq;
    Mq(:,:,2,k) = Yq;
    Mq(:,:,3,k) = Zq;
    
end




end

