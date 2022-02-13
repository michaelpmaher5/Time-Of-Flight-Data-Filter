function RawValAve = rawvalueave(Xraw,Yraw,Zraw,nbrOfFrames)

    RawValAve = zeros(size(Xraw));
    
    for x = 1:size(Xraw,1)
        for y = 1:size(Xraw,2)
            for k = 1:nbrOfFrames-2
                RawValAve(x,y,k) = sqrt((Xraw(x,y,k)^2)+(Yraw(x,y,k)^2)+(Zraw(x,y,k)^2));
            end
        end
    end
end