function avestddev = getAverageStandardDeviation(stddev)
    avestddev = 0;
    
    for x = 1:size(stddev,1)
        for y = 1:size(stddev,2)
            avestddev = avestddev + stddev(x,y,3);
        end
    end
    avestddev = avestddev/(size(stddev,1) * size(stddev,2));
end