function [Vraw,Araw] = raw_time_derivatives(rawpos,fps,nbrOfFrames)
    Vraw = zeros(size(rawpos,1),size(rawpos,2),nbrOfFrames-1);
    Araw = zeros(size(rawpos,1),size(rawpos,2),nbrOfFrames-2);

    h = 1/fps;    %% Original value: 1/2*1/fps (delta t)
    
    for t=1:nbrOfFrames-1
    
    Vraw(:,:,t) = (rawpos(:,:,t+1)-rawpos(:,:,t))/h;
   
    
    end

for t=1:nbrOfFrames-2
    %{
    A(:,:,1,t) = (M(:,:,1,t+2)+M(:,:,1,t+1)-2*M(:,:,1,t))/h^2;
    A(:,:,2,t) = (M(:,:,2,t+2)+M(:,:,2,t+1)-2*M(:,:,2,t))/h^2;
    A(:,:,3,t) = (M(:,:,3,t+2)+M(:,:,3,t+1)-2*M(:,:,3,t))/h^2;
   %}
    
    Araw(:,:,t) = (Vraw(:,:,t+1)-Vraw(:,:,t))/h;
   
end
end