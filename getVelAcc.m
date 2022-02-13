function V = getVelAcc(filtered,fps,nbrOfFrames)
    V = zeros(size(filtered));
    t = 1/fps;    %% Original value: 1/2*1/fps (delta t)
    
    %% Initial Conditions
    for x = 1:size(filtered,1)
        for y = 1:size(filtered,2)
            V(x,y,1,1) = 0;           
            V(x,y,2,1) = 0;           
            V(x,y,3,1) = 0;           
        end
    end
    
    %% Calculate Velocity and/or Acceleration
   for x = 1:size(filtered,1)
       for y = 1:size(filtered,2)
           for k = 2:nbrOfFrames
               V(x,y,1,k) = (filtered(x,y,1,k) - filtered(x,y,1,k-1))/t;
               V(x,y,2,k) = (filtered(x,y,2,k) - filtered(x,y,2,k-1))/t;
               V(x,y,3,k) = (filtered(x,y,3,k) - filtered(x,y,3,k-1))/t;
           end
       end
   end
end