filepath = "Earth_nadir/";
refine_factor = 4;
load(filepath + "Ampraw.mat")

%Ampraw = Ampraw_avg;

[pxX,pxY,nbrOfFrames] = size(Ampraw);
%Vq = zeros((pxX-1)*2^refine_factor+1,(pxY-1)*2^refine_factor+1,1);



%%% Frame average
nbrOfFrames2 = floor(nbrOfFrames / 10) + min(mod(nbrOfFrames,10),1);
Ampraw_avg = zeros(pxX,pxY,nbrOfFrames2);

for k=1:nbrOfFrames2
    
    Beg = (k-1)*10+1;
    End = min(k*10,nbrOfFrames);
    
    Ampraw_avg(:,:,k) = mean_global_matrix(Ampraw(:,:,Beg:End),numel(Beg:End));
    
end



Vq(:,:) = interp2(Ampraw(:,:,1),refine_factor,'makima');
Vq_avg(:,:) = interp2(Ampraw_avg(:,:,1),refine_factor,'makima');


%%% Re-scale
Ampraw(:,:,1) = Ampraw(:,:,1) - min(min(min(Ampraw(:,:,1))));
Ampraw(:,:,1) = Ampraw(:,:,1)/(max(max(max(Ampraw(:,:,1)))));
Vq = Vq - min(min(min(Vq)));
Vq = Vq/(max(max(max(Vq))));
Ampraw_avg(:,:,1) = Ampraw_avg(:,:,1) - min(min(min(Ampraw_avg(:,:,1))));
Ampraw_avg(:,:,1) = Ampraw_avg(:,:,1)/(max(max(max(Ampraw_avg(:,:,1)))));
Vq_avg = Vq_avg - min(min(min(Vq_avg)));
Vq_avg = Vq_avg/(max(max(max(Vq_avg))));


%%% Images
imwrite(Ampraw(:,:,1),filepath+'single_image.png');
imwrite(Vq,filepath+'single_image_interpolated.png');
imwrite(Ampraw_avg(:,:,1),filepath+'all_frames_averaged_image.png');
imwrite(Vq_avg,filepath+'all_frames_averaged_image_interpolated.png');



%%% Figures
figure('Units','normalized','Position',[0 0 1 1]); 
imagesc(Ampraw(:,:,1));
colormap gray
axis image
axis off
title('Original Image');

figure('Units','normalized','Position',[0 0 1 1]); 
imagesc(Vq(:,:,1));
colormap gray
axis image
axis off
title('Interpolation');

figure('Units','normalized','Position',[0 0 1 1]); 
imagesc(Ampraw_avg(:,:,1));
colormap gray
axis image
axis off
title('Averaged Image');

figure('Units','normalized','Position',[0 0 1 1]); 
imagesc(Vq_avg(:,:,1));
colormap gray
axis image
axis off
title('Averaged Interpolation');