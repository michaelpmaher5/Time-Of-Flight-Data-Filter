function [ Std ] = std_matrix( M, nbrOfFrames )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if numel(size(M))==3

    Std = zeros(size(M,1),size(M,2));

    for i=1:size(M,1)
        for j=1:size(M,2)

            Xtemp=[];

            for k=1:nbrOfFrames
                if (~isnan(M(i,j,k)))
                    Xtemp=[Xtemp M(i,j,k)];
                end
            end

            Std(i,j)=std(Xtemp,1);
            
        end
    end
end

if numel(size(M))==4

    Std = zeros(size(M,1),size(M,2),3);

    for i=1:size(M,1)
        for j=1:size(M,2)

            Xtemp=[];
            Ytemp=[];
            Ztemp=[];

            for k=1:nbrOfFrames

                if (~isnan(M(i,j,1,k)))
                    Xtemp=[Xtemp M(i,j,1,k)];
                end
                if (~isnan(M(i,j,2,k)))
                    Ytemp=[Ytemp M(i,j,2,k)];
                end
                if (~isnan(M(i,j,3,k)))
                    Ztemp=[Ztemp M(i,j,3,k)];
                end

            end

            Std(i,j,1)=std(Xtemp,1);
            Std(i,j,2)=std(Ytemp,1);
            Std(i,j,3)=std(Ztemp,1);
        end
    end

end

end

