function normalsoriginal = getNormalCoordinates(sail,nbrOfFrames)


normalsoriginal = zeros(size(sail));


for k = 1:nbrOfFrames
    [Nx,Ny,Nz] = surfnorm(sail(:,:,1,k),sail(:,:,2,k),sail(:,:,3,k));
    for x = 1:size(sail,1)
        for y = 1:size(sail,2)
            Nx(x,y) = Nx(x,y) - sail(x,y,1,k);
            Ny(x,y) = Ny(x,y) - sail(x,y,2,k);
            Nz(x,y) = Nz(x,y) - sail(x,y,3,k);
        end
    end
    for x1 = 1:size(Nx,1)
        for y1 = 1:size(Nx,2)
            normalsoriginal(x1,y1,1,k) = Nx(x1,y1);
        end
    end
    
     for x1 = 1:size(Ny,1)
        for y1 = 1:size(Ny,2)
            normalsoriginal(x1,y1,2,k) = Ny(x1,y1);
        end
     end
    
      for x1 = 1:size(Nz,1)
        for y1 = 1:size(Nz,2)
            normalsoriginal(x1,y1,3,k) = Nz(x1,y1);
        end
      end
end 

   
   
end