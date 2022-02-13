function [ Mean ] = mean_global_matrix( M, nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if numel(size(M))==3

    Mean = zeros(size(M,1),size(M,2));

    for i=1:size(M,1)
        for j=1:size(M,2)

            Xtemp=[];

            for k=1:nbrOfFrames
                if (~isnan(M(i,j,k)))
                    Xtemp=[Xtemp M(i,j,k)];
                end
            end

            Mean(i,j)=mean(Xtemp);
            
        end
    end
end

if numel(size(M))==4

    Mean = zeros(size(M,1),size(M,2),3);

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

            Mean(i,j,1)=mean(Xtemp);
            Mean(i,j,2)=mean(Ytemp);
            Mean(i,j,3)=mean(Ztemp);
        end
    end
end



end

