function [ a,b,c,d ] = plane( S1,S2,S3,S4,nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = zeros(nbrOfFrames,1);
b = zeros(nbrOfFrames,1);
c = zeros(nbrOfFrames,1);
d = zeros(nbrOfFrames,1);

for k=1:nbrOfFrames
    
    n = mean([ cross(S2(:,1,k)-S1(:,1,k),S4(:,1,k)-S1(:,1,k)) , cross(S3(:,1,k)-S2(:,1,k),S1(:,1,k)-S2(:,1,k)) , cross(S4(:,1,k)-S3(:,1,k),S2(:,1,k)-S3(:,1,k)) , cross(S1(:,1,k)-S4(:,1,k),S3(:,1,k)-S4(:,1,k)) ] , 2);
    
    a(k) = n(1);
    b(k) = n(2);
    c(k) = n(3);
    d(k) = -( a(k)*mean([S1(1,1,k),S2(1,1,k),S3(1,1,k),S4(1,1,k)]) + b(k)*mean([S1(2,1,k),S2(2,1,k),S3(2,1,k),S4(2,1,k)]) + c(k)*mean([S1(3,1,k),S2(3,1,k),S3(3,1,k),S4(3,1,k)]) );
    
end

end

