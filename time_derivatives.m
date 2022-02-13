function [ V,A ] = time_derivatives( M,fps,nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


V = zeros(size(M,1),size(M,2),3,nbrOfFrames-1);
A = zeros(size(M,1),size(M,2),3,nbrOfFrames-2);
%%M = M/1000;
%%h = 1/2*1/fps;
h = 1/fps;    %% Original value: 1/2*1/fps (delta t)

for t=1:nbrOfFrames-1
    
    V(:,:,1,t) = (M(:,:,1,t+1)-M(:,:,1,t))/h;
    V(:,:,2,t) = (M(:,:,2,t+1)-M(:,:,2,t))/h;
    V(:,:,3,t) = (M(:,:,3,t+1)-M(:,:,3,t))/h;
    
end

for t=1:nbrOfFrames-2
    %{
    A(:,:,1,t) = (M(:,:,1,t+2)+M(:,:,1,t+1)-2*M(:,:,1,t))/h^2;
    A(:,:,2,t) = (M(:,:,2,t+2)+M(:,:,2,t+1)-2*M(:,:,2,t))/h^2;
    A(:,:,3,t) = (M(:,:,3,t+2)+M(:,:,3,t+1)-2*M(:,:,3,t))/h^2;
   %}
    
    A(:,:,1,t) = (V(:,:,1,t+1)-V(:,:,1,t))/h;
    A(:,:,2,t) = (V(:,:,2,t+1)-V(:,:,2,t))/h;
    A(:,:,3,t) = (V(:,:,3,t+1)-V(:,:,3,t))/h;
end

end
