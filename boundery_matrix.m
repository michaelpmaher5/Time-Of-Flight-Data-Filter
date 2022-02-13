function [ M_bound, A1, A2, A3, A4 ] = boundery_matrix( S1, S2, S3, S4, nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


M_bound = zeros(14,14,2,nbrOfFrames);

A1=zeros(2,14,nbrOfFrames);
A2=zeros(2,14,nbrOfFrames);
A3=zeros(2,14,nbrOfFrames);
A4=zeros(2,14,nbrOfFrames);

vectX=zeros(2,14,nbrOfFrames);
vectY=zeros(2,14,nbrOfFrames);

for k=1:nbrOfFrames
    
    A1(:,1,k) = [S1(1,1,k)-(S2(1,1,k)-S1(1,1,k))/24-(S3(1,1,k)-S2(1,1,k))/24 ; S1(2,1,k)-(S2(2,1,k)-S1(2,1,k))/24-(S3(2,1,k)-S2(2,1,k))/24];
    A2(:,1,k) = [S2(1,1,k)-(S3(1,1,k)-S2(1,1,k))/24+(S2(1,1,k)-S1(1,1,k))/24 ; S2(2,1,k)-(S3(2,1,k)-S2(2,1,k))/24+(S2(2,1,k)-S1(2,1,k))/24];
    A3(:,1,k) = [S4(1,1,k)-(S3(1,1,k)-S4(1,1,k))/24+(S4(1,1,k)-S1(1,1,k))/24 ; S4(2,1,k)-(S3(2,1,k)-S4(2,1,k))/24+(S4(2,1,k)-S1(2,1,k))/24];
    A4(:,1,k) = [S1(1,1,k)-(S4(1,1,k)-S1(1,1,k))/24-(S3(1,1,k)-S4(1,1,k))/24 ; S1(2,1,k)-(S4(2,1,k)-S1(2,1,k))/24-(S3(2,1,k)-S4(2,1,k))/24];
    
    for i=2:14
        A1(:,i,k) = [A1(1,i-1,k)+(S2(1,1,k)-S1(1,1,k))/12 ; A1(2,i-1,k)+(S2(2,1,k)-S1(2,1,k))/12];
        A2(:,i,k) = [A2(1,i-1,k)+(S3(1,1,k)-S2(1,1,k))/12 ; A2(2,i-1,k)+(S3(2,1,k)-S2(2,1,k))/12];
        A3(:,i,k) = [A3(1,i-1,k)+(S3(1,1,k)-S4(1,1,k))/12 ; A3(2,i-1,k)+(S3(2,1,k)-S4(2,1,k))/12];
        A4(:,i,k) = [A4(1,i-1,k)+(S4(1,1,k)-S1(1,1,k))/12 ; A4(2,i-1,k)+(S4(2,1,k)-S1(2,1,k))/12];
    end

    
    for p=1:14
        vectX(:,p,k) = [(A3(1,p,k)-A1(1,p,k))/13 ; (A3(2,p,k)-A1(2,p,k))/13];
        vectY(:,p,k) = [(A2(1,p,k)-A4(1,p,k))/13 ; (A2(2,p,k)-A4(2,p,k))/13];
    end
    
    
    for i=1:14
        for j=1:14
            
            M_bound(i,j,1,k) = A1(1,i,k) + (j-1)*vectX(1,i,k) ;     %+ (i-1)*vectY(1,j,k);
            M_bound(i,j,2,k) = A4(2,j,k) + (i-1)*vectY(2,j,k) ;     %+ (j-1)*vectX(2,i,k);
            
        end
    end
    
end



end

