function [Xrawposmean,Yrawposmean,Zrawposmean] = rawinterpolate(Xraw,Yraw,Zraw,nbrOfFrames)

Xrawposmean = zeros(size(Xraw));
Yrawposmean = zeros(size(Yraw));
Zrawposmean = zeros(size(Zraw));

coarse = 5;

%% Fill Xrawposmean
for x = 1:size(Xraw,1)
    for y = 1:size(Xraw,2)
        for k = 1:coarse
            Xrawposmean(x,y,k) = (Xraw(x,y,k)+Xraw(x,y,k+1)+Xraw(x,y,k+2)+Xraw(x,y,k+3)+Xraw(x,y,k+4))/coarse;
        end
        for k2 = nbrOfFrames-coarse:nbrOfFrames
            Xrawposmean(x,y,k2) = (Xraw(x,y,k2)+Xraw(x,y,k2-1)+Xraw(x,y,k2-2)+Xraw(x,y,k2-3)+Xraw(x,y,k2-4))/coarse;
        end
        for k3 = coarse+1:nbrOfFrames-coarse-1
            Xrawposmean(x,y,k3) = (Xraw(x,y,k3-2)+Xraw(x,y,k3-1)+Xraw(x,y,k3)+Xraw(x,y,k3+1)+Xraw(x,y,k3+2))/coarse;
        end
    end
end

%% Fill Yrawposmean
for x = 1:size(Yraw,1)
    for y = 1:size(Yraw,2)
        for k = 1:coarse
            Yrawposmean(x,y,k) = (Yraw(x,y,k)+Yraw(x,y,k+1)+Yraw(x,y,k+2)+Yraw(x,y,k+3)+Yraw(x,y,k+4))/coarse;
        end
        for k2 = nbrOfFrames-coarse:nbrOfFrames
            Yrawposmean(x,y,k2) = (Yraw(x,y,k2)+Yraw(x,y,k2-1)+Yraw(x,y,k2-2)+Yraw(x,y,k2-3)+Yraw(x,y,k2-4))/coarse;
        end
        for k3 = coarse+1:nbrOfFrames-coarse-1
            Yrawposmean(x,y,k3) = (Yraw(x,y,k3-2)+Yraw(x,y,k3-1)+Yraw(x,y,k3)+Yraw(x,y,k3+1)+Yraw(x,y,k3+2))/coarse;
        end
    end
end

%% Fill Zrawposmean
for x = 1:size(Zraw,1)
    for y = 1:size(Zraw,2)
        for k = 1:coarse
            Zrawposmean(x,y,k) = (Zraw(x,y,k)+Zraw(x,y,k+1)+Zraw(x,y,k+2)+Zraw(x,y,k+3)+Zraw(x,y,k+4))/coarse;
        end
        for k2 = nbrOfFrames-coarse:nbrOfFrames
            Zrawposmean(x,y,k2) = (Zraw(x,y,k2)+Zraw(x,y,k2-1)+Zraw(x,y,k2-2)+Zraw(x,y,k2-3)+Zraw(x,y,k2-4))/coarse;
        end
        for k3 = coarse+1:nbrOfFrames-coarse-1
            Zrawposmean(x,y,k3) = (Zraw(x,y,k3-2)+Zraw(x,y,k3-1)+Zraw(x,y,k3)+Zraw(x,y,k3+1)+Zraw(x,y,k3+2))/coarse;
        end
    end
end
end
