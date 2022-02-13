function [ S1,S2,S3,S4,L1,L2,L3,L4 ] = summit( Z, nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% We get the outer corners
S1=zeros(3,1,nbrOfFrames);
S2=zeros(3,1,nbrOfFrames);
S3=zeros(3,1,nbrOfFrames);
S4=zeros(3,1,nbrOfFrames);

s11=zeros(3,1,nbrOfFrames);
s12=zeros(3,1,nbrOfFrames);
s21=zeros(3,1,nbrOfFrames);
s22=zeros(3,1,nbrOfFrames);
s31=zeros(3,1,nbrOfFrames);
s32=zeros(3,1,nbrOfFrames);
s41=zeros(3,1,nbrOfFrames);
s42=zeros(3,1,nbrOfFrames);

L1=zeros(nbrOfFrames,1);
L2=zeros(nbrOfFrames,1);
L3=zeros(nbrOfFrames,1);
L4=zeros(nbrOfFrames,1);


for k=1:nbrOfFrames
    
a1=2;   a2=2;   a3=2;    a4=2;
    
    %S1
    for i=1:1:30
        for j=30:1:60
            if Z(i,j,k)~=0
                s11(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s11(:,:,k)~=[0;0;0]
            break;
        end
    end
    for j=30:1:60
        for i=1:1:30
            if Z(i,j,k)~=0
                s12(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s12(:,:,k)~=[0;0;0]
            break;
        end
    end
    S1(:,:,k) = 1/2*(s11(:,:,k)+s12(:,:,k));
    S1(1,:,k) = floor(S1(1,:,k));
    S1(2,:,k) = floor(S1(2,:,k));
    
    while a1>0
        a1=2;
        if Z(S1(2,:,k),S1(1,:,k)-1,k)==0
            a1=a1-1;
        else
            S1(1,:,k)=S1(1,:,k)-1;
        end
        if Z(S1(2,:,k)-1,S1(1,:,k),k)==0
            a1=a1-1;
        else
            S1(2,:,k)=S1(2,:,k)-1;
        end
    end    
    
    
    %S2
    for i=110:-1:80
        for j=30:1:60
            if Z(i,j,k)~=0
                s21(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s21(:,:,k)~=[0;0;0]
            break;
        end
    end
    for j=30:1:60
        for i=110:-1:80
            if Z(i,j,k)~=0
                s22(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s22(:,:,k)~=[0;0;0]
            break;
        end
    end
    S2(:,:,k) = 1/2*(s21(:,:,k)+s22(:,:,k));
    S2(1,:,k) = floor(S2(1,:,k));
    S2(2,:,k) = ceil(S2(2,:,k));
    
    while a2>0
        a2=2;
        if Z(floor(S2(2,:,k)),S2(1,:,k)-1,k)==0
            a2=a2-1;
        else
            S2(1,:,k)=S2(1,:,k)-1;
        end
        if Z(S2(2,:,k)+1,S2(1,:,k),k)==0
            a2=a2-1;
        else
            S2(2,:,k)=S2(2,:,k)+1;
        end
    end    
    
    
    %S3
    for i=110:-1:80
        for j=140:-1:110
            if Z(i,j,k)~=0
                s31(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s31(:,:,k)~=[0;0;0]
            break;
        end
    end
    for j=140:-1:110
        for i=110:-1:80
            if Z(i,j,k)~=0
                s32(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s32(:,:,k)~=[0;0;0]
            break;
        end
    end
    S3(:,:,k) = 1/2*(s31(:,:,k)+s32(:,:,k));
    S3(1,:,k) = ceil(S3(1,:,k));
    S3(2,:,k) = ceil(S3(2,:,k));
    
    while a3>0
        a3=2;
        if Z(S3(2,:,k),S3(1,:,k)+1,k)==0
            a3=a3-1;
        else
            S3(1,:,k)=S3(1,:,k)+1;
        end
        if Z(S3(2,:,k)+1,S3(1,:,k),k)==0
            a3=a3-1;
        else
            S3(2,:,k)=S3(2,:,k)+1;
        end
    end    
    
    
    %S4
    for i=1:1:30
        for j=140:-1:110
            if Z(i,j,k)~=0
                s41(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s41(:,:,k)~=[0;0;0]
            break;
        end
    end
    for j=140:-1:110
        for i=1:1:30
            if Z(i,j,k)~=0
                s42(:,:,k) = [j;i;Z(i,j,k)];
                break;
            end
        end
        if s42(:,:,k)~=[0;0;0]
            break;
        end
    end
    S4(:,:,k) = 1/2*(s41(:,:,k)+s42(:,:,k));
    S4(1,:,k) = ceil(S4(1,:,k));
    S4(2,:,k) = floor(S4(2,:,k));
    
    while a4>0
        a4=2;
        if Z(S4(2,:,k),S4(1,:,k)+1,k)==0
            a4=a4-1;
        else
            S4(1,:,k)=S4(1,:,k)+1;
        end
        if Z(S4(2,:,k)-1,S4(1,:,k),k)==0
            a4=a4-1;
        else
            S4(2,:,k)=S4(2,:,k)-1;
        end
    end    
    
    
    
    
    
    
    
    %L_i
    L1(k) = sqrt((S2(1,1,k)-S1(1,1,k))^2 + (S2(2,1,k)-S1(2,1,k))^2);
    L2(k) = sqrt((S3(1,1,k)-S2(1,1,k))^2 + (S3(2,1,k)-S2(2,1,k))^2);
    L3(k) = sqrt((S4(1,1,k)-S3(1,1,k))^2 + (S4(2,1,k)-S3(2,1,k))^2);
    L4(k) = sqrt((S1(1,1,k)-S4(1,1,k))^2 + (S1(2,1,k)-S4(2,1,k))^2);

end

end

