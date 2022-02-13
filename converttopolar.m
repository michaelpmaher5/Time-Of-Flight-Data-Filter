function polar = converttopolar(cart,nbrOfFrames)
    polar = zeros(size(cart));
    
    for x = 1:size(cart,1)
        for y = 1:size(cart,2)
            for k = 1:nbrOfFrames
                polar(x,y,1,k) = sqrt(cart(x,y,1,k)^2 + cart(x,y,2,k)^2 + cart(x,y,3,k)^2);
                polar(x,y,2,k) = acos(cart(x,y,3,k)/polar(x,y,1,k));
                polar(x,y,3,k) = acos(cart(x,y,1,k)/(polar(x,y,1,k) * sin(polar(x,y,2,k))));
            end
        end
    end
    
end
