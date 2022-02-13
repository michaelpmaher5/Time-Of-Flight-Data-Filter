function stddev = standarddeviation(sailmatrix,nbrOfFrames)
    
    stddev = zeros(size(sailmatrix,1),size(sailmatrix,2),3);
    mean = zeros(size(sailmatrix,1),size(sailmatrix,2),3);
    
    %% First, find the mean z values at each pixel and store them in the mean matrix
    for x = 1:size(sailmatrix,1)
        for y = 1:size(sailmatrix,2)
            for k = 1:nbrOfFrames
                mean(x,y,1) = mean(x,y,1) + sailmatrix(x,y,1,k);
                mean(x,y,2) = mean(x,y,2) + sailmatrix(x,y,2,k);
                mean(x,y,3) = mean(x,y,3) + sailmatrix(x,y,3,k);
            end
            mean(x,y,1) = mean(x,y,1)/nbrOfFrames;
            mean(x,y,2) = mean(x,y,2)/nbrOfFrames;
            mean(x,y,3) = mean(x,y,3)/nbrOfFrames;
        end
    end

    %% Next, calculate the standard deviation at each pixel over all the frames
    for x = 1:size(sailmatrix,1)
        for y = 1:size(sailmatrix,2)
            for k = 1:nbrOfFrames
                stddev(x,y,1) = stddev(x,y,1) + (sailmatrix(x,y,1,k) - mean(x,y,1))^2;
                stddev(x,y,2) = stddev(x,y,2) + (sailmatrix(x,y,2,k) - mean(x,y,2))^2;
                stddev(x,y,3) = stddev(x,y,3) + (sailmatrix(x,y,3,k) - mean(x,y,3))^2;
            end
            stddev(x,y,1) = (stddev(x,y,1)/(nbrOfFrames-1))^0.5;
            stddev(x,y,2) = (stddev(x,y,2)/(nbrOfFrames-1))^0.5;
            stddev(x,y,3) = (stddev(x,y,3)/(nbrOfFrames-1))^0.5;
        end
    end
end
    
