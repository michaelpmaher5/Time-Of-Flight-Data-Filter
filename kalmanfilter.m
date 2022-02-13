function filtered = kalmanfilter(sailmatrix,stddev,nbrOfFrames)

    filtered = zeros(size(sailmatrix));
    estimateuncertainty = zeros(size(sailmatrix));
    kn = zeros(size(sailmatrix));         %% kalman gains
    
    %% Initial Conditions (assuming initial values are correct in matrix)
    for x = 1:size(sailmatrix,1)
        for y = 1:size(sailmatrix,2)
            estimateuncertainty(x,y,1,1) = stddev(x,y,1);
            filtered(x,y,1,1) = sailmatrix(x,y,1,1);
            estimateuncertainty(x,y,2,1) = stddev(x,y,2);
            filtered(x,y,2,1) = sailmatrix(x,y,2,1);
            estimateuncertainty(x,y,3,1) = stddev(x,y,3);
            filtered(x,y,3,1) = sailmatrix(x,y,3,1);
        end
    end
    
    %% Filter Out the Noise
    for x = 1:size(sailmatrix,1)
        for y = 1:size(sailmatrix,2)
            for k = 2:nbrOfFrames
                % X Values
                kn(x,y,1,k) = estimateuncertainty(x,y,1,k-1)/(estimateuncertainty(x,y,1,k-1) + stddev(x,y,1));
                estimateuncertainty(x,y,1,k) = (1 - kn(x,y,1,k))*estimateuncertainty(x,y,1,k-1);
                filtered(x,y,1,k) = filtered(x,y,1,k-1) + (kn(x,y,1,k))*(sailmatrix(x,y,1,k) - filtered(x,y,1,k-1));
                % Y Values
                kn(x,y,2,k) = estimateuncertainty(x,y,2,k-1)/(estimateuncertainty(x,y,2,k-1) + stddev(x,y,2));
                estimateuncertainty(x,y,2,k) = (1 - kn(x,y,2,k))*estimateuncertainty(x,y,2,k-1);
                filtered(x,y,2,k) = filtered(x,y,2,k-1) + (kn(x,y,2,k))*(sailmatrix(x,y,2,k) - filtered(x,y,2,k-1));
                %Z Values
                kn(x,y,3,k) = estimateuncertainty(x,y,3,k-1)/(estimateuncertainty(x,y,3,k-1) + stddev(x,y,3));
                estimateuncertainty(x,y,3,k) = (1 - kn(x,y,3,k))*estimateuncertainty(x,y,3,k-1);
                filtered(x,y,3,k) = filtered(x,y,3,k-1) + (kn(x,y,3,k))*(sailmatrix(x,y,3,k) - filtered(x,y,3,k-1));
            end   
        end
    end
    
end